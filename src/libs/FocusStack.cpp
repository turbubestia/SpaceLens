/*
 * FocusStack.cpp
 *
 *  Created on: Aug 29, 2020
 *      Author: claud
 */

#include <QtCore/QtCore>

#include <opencv4/opencv2/opencv.hpp>
#include <opencv4/opencv2/xfeatures2d.hpp>
#include <opencv4/opencv2/ximgproc/edge_filter.hpp>

#include <eigen3/Eigen/Eigen>

#include <Utils.h>

#include "FocusStack.h"

cv::Mat vec2f2map(cv::Mat x, cv::Mat y) {

    double minX, minY;
    double maxX, maxY;
    minMaxLoc( x, &minX, &maxX, nullptr, nullptr );
    minMaxLoc( y, &minY, &maxY, nullptr, nullptr );
    double maxVal = fmax(fabs(maxX), fabs(maxY));
    double minVal = fmin(fabs(minX), fabs(minY));
    // range defined as [-range,0,+range] interval
    double range = 2.0*fmax(maxVal, minVal);
    // Map [-range,0,+range] to [0,1]
    cv::Mat x_temp = x/range + 0.5;
    cv::Mat y_temp = y/range + 0.5;
    // Map an float image (x,y) to a BGR image (b,0,r)
    cv::Mat colormap(x.size(), CV_8UC3);

    for (int j = 0; j < x.rows; j++)
    {
        float* x_ptr = x_temp.ptr<float>(j);
        float* y_ptr = y_temp.ptr<float>(j);

        cv::Vec3b *colormap_ptr = colormap.ptr<cv::Vec3b>(j);

        for (int k = 0; k < x.cols; k++) {
            colormap_ptr[k][0] = y_ptr[k]*255; // B
            colormap_ptr[k][1] = 0; // G
            colormap_ptr[k][2] = x_ptr[k]*255; // R
        }
    }

    return colormap;
}

cv::Mat depth2map(cv::Mat img) {
    // range defined as [0,maxVal]
    double maxVal;
    minMaxLoc( img, nullptr, &maxVal, nullptr, nullptr );
    // Map an integer image to BGR rainbow image
    cv::Mat colormap(img.size(), CV_8UC3);

    for (int j = 0; j < img.rows; j++) {
        uchar* temp_ptr = img.ptr<uchar>(j);
        cv::Vec3b *colormap_ptr = colormap.ptr<cv::Vec3b>(j);

        for (int k = 0; k < img.cols; k++) {
            colormap_ptr[k][0] = temp_ptr[k] * 150 / maxVal;
            colormap_ptr[k][1] = 227;
            colormap_ptr[k][2] = temp_ptr[k] == 0 ? 0 : 227;
        }
    }

    // convert from Hue to BGR
    cvtColor( colormap, colormap, cv::COLOR_HSV2BGR);

    return colormap;
}

cv::Mat scalar2gray(cv::Mat img) {
    // range defined as [0,maxVal]
    double maxVal;
    minMaxLoc( img, nullptr, &maxVal, nullptr, nullptr );
    // Map a float image to BGR grayscale image
    cv::Mat colormap(img.size(), CV_8UC1);

    for (int j = 0; j < img.rows; j++) {
        float* gray_ptr = img.ptr<float>(j);
        uchar* colormap_ptr = colormap.ptr<uchar>(j);

        for (int k = 0; k < img.cols; k++) {
            colormap_ptr[k] = gray_ptr[k] * 255 / maxVal;
        }
    }

    cvtColor( colormap, colormap, cv::COLOR_GRAY2BGR);

    return colormap;
}

void downSampleToGray(cv::Mat &src, cv::Mat &dst) {
    // downsample by a factor of 2
    dst = cv::Mat(src.size()/2, CV_8UC1);
    // gaussian kernel
    float data[3][3] = {{1,2,1},{2,4,2},{1,2,1}};

    for (int i = 1; i < dst.rows-1; i++)
    {
        // fast access to y-1, y, and y+1 pixels
        cv::Vec3b *src_ptr[3];
        src_ptr[0] = src.ptr<cv::Vec3b>(2*i-1);
        src_ptr[1] = src.ptr<cv::Vec3b>(2*i);
        src_ptr[2] = src.ptr<cv::Vec3b>(2*i+1);

        uchar *dst_ptr = dst.ptr<uchar>(i);

        for (int j = 1; j < dst.cols-1; j++)
        {
            // In one step compute the gaussian average and the equivalent
            // grayscale of the downsample image
            float avg = 0;
            for(int m = 0; m < 3; m++) {
                for(int n = -1; n <= 1; n++) {
                    // pixel color
                    cv::Vec3b bgr = src_ptr[m][2*j+n];
                    // average contribution
                    avg += (0.299*bgr[2] + 0.587*bgr[1] + 0.114*bgr[0]) * data[m][n];
                }
            }

            // Normalization of the gaussian kernel
            dst_ptr[j] = avg/16.0;
        }
    }
}

FocusStack::FocusStack(QObject *parent)
    : QObject(parent) {
    
    int nbWorkers = QThread::idealThreadCount();

    for (int i = 0; i < nbWorkers; i++) {
        FocusStackWorker *worker = new FocusStackWorker;
        connect(worker, &FocusStackWorker::imageFocused, this, &FocusStack::imageFocused);
        connect(worker, &FocusStackWorker::imageAligned, this, &FocusStack::imageAligned);
        connect(worker, &FocusStackWorker::workerFinished, this, &FocusStack::workerFinished);
        workerList.append(worker);
    }

    finishedJobs = 0;
    workersSpanwned = 0;
}

FocusStack::~FocusStack() {

}

void FocusStack::addStackGroup(QStringList stack) {
    groups.append(stack);
}

void FocusStack::clear() {
    groups.clear();

    for (int i = 0; i < workerList.count(); i++) {
        workerList.at(i)->clear();
    }

    finishedJobs = 0;
    workersSpanwned = 0;
}

void FocusStack::setCacheFolder(const QString &folder) {
    cacheFolder = folder;
}

void FocusStack::setOutputFolder(const QString &folder) {
    outputFolder = folder;
}

void FocusStack::startFocusStack() {
    // divide all the groups to spread and fit into the workerlist
    int nbWorkers = workerList.count();
    int groupSize = ceil((double) groups.count() / nbWorkers);
    if (groupSize < 1) {
        groupSize = 1;
    }

    finishedJobs = 0;
    workersSpanwned = 0;

    for (int w = 0, g = 0; w < nbWorkers && g < groups.count(); w++)
    {
        // adjust the remaining groups
        int currentGroupSize = groupSize;
        int remElements = groups.count() - g;
        if (currentGroupSize > remElements) {
            currentGroupSize = remElements;
        }
        // setup the worker
        workerList.at(w)->setCacheFolder(cacheFolder);
        workerList.at(w)->setOutputFolder(outputFolder);
        for(int j = 0; j < currentGroupSize; j++, g++) {
            workerList.at(w)->addGroup(g, groups.at(g));
        }
        workerList.at(w)->start();
        // keep track of the groups already spawned
        workersSpanwned++;
    }
}

void FocusStack::workerFinished() {
    finishedJobs++;
    if (finishedJobs == workersSpanwned) {
        emit finished();
    }
}

FocusStackWorker::FocusStackWorker(QObject *parent)
    : QThread(parent) {

    currentGroup = 0;
}

void FocusStackWorker::align() {

    QStringList filenames = groups.at(currentGroup);

    // The first images is the reference
    cv::Mat ref_color = cv::imread(filenames.at(0).toStdString());
    // it is aligned to itself
    cv::imwrite(cacheFilenames.at(0).toStdString(), ref_color);
    // use the downsampled image as reference for fast aligning
    cv::Mat ref_gray;
    downSampleToGray(ref_color, ref_gray);

    emit imageAligned(groupIds.at(currentGroup), filenames.at(0));

    // Align the rest of images to the reference
    for (int i = 1; i < filenames.size(); i++) {
        // Set the warp matrix to identity.
        cv::Mat warp_matrix = cv::Mat::eye(2, 3, CV_32F);

        // Set the stopping criteria for the algorithm.
        int number_of_iterations = 200;
        double termination_eps = 1e-6;
        cv::TermCriteria criteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, number_of_iterations, termination_eps);

        // next image to align to the reference
        cv::Mat src_color;
        src_color = cv::imread(filenames.at(i).toStdString());

        // Find the transformation using the downsampled version of the image
        cv::Mat src_gray;
        downSampleToGray(src_color, src_gray);
        cv::findTransformECC(ref_gray, src_gray, warp_matrix, cv::MOTION_AFFINE, criteria);

        // The current image is the reference for the next one
        cv::warpAffine(src_gray, ref_gray, warp_matrix, src_gray.size(), cv::INTER_LINEAR + cv::WARP_INVERSE_MAP);

        // Scale the non-linear factors of the affine matrix to account for the half-size of
        // the reference image used to find the transformation
        warp_matrix.at<float>(0,2) *= 2;
        warp_matrix.at<float>(1,2) *= 2;

        // Align the source image
        cv::Mat src_aligned;
        cv::warpAffine(src_color, src_aligned, warp_matrix, src_color.size(), cv::INTER_LINEAR + cv::WARP_INVERSE_MAP);

        // Save aligned
        cv::imwrite(cacheFilenames.at(i).toStdString(), src_aligned);

        emit imageAligned(groupIds.at(currentGroup), filenames.at(i));
    }
}

void colorDerivative(int depth, cv::Mat src, cv::Mat *dx, cv::Mat *dy, cv::Mat &M, cv::Mat &D) {
    cv::Mat img, channels[3];

    // Low pass filter
    cv::GaussianBlur(src, img, cv::Size(3,3), -1);
    // Separate the B,G, and R channels
    cv::split(img, channels);

    // Get the derivative of each color to get the most information out of it, as
    // two different colors can have the same gray equivalence
    cv::Scharr(channels[0], dx[0], CV_32F, 1, 0);
    cv::Scharr(channels[1], dx[1], CV_32F, 1, 0);
    cv::Scharr(channels[2], dx[2], CV_32F, 1, 0);

    cv::Scharr(channels[0], dy[0], CV_32F, 0, 1);
    cv::Scharr(channels[1], dy[1], CV_32F, 0, 1);
    cv::Scharr(channels[2], dy[2], CV_32F, 0, 1);

    for (int j = 0; j < M.rows; j++)
    {
        // fast access pointers
        float* dx0_ptr = dx[0].ptr<float>(j);
        float* dx1_ptr = dx[1].ptr<float>(j);
        float* dx2_ptr = dx[2].ptr<float>(j);

        float* dy0_ptr = dy[0].ptr<float>(j);
        float* dy1_ptr = dy[1].ptr<float>(j);
        float* dy2_ptr = dy[2].ptr<float>(j);

        // M holds the maximum gradient
        float* M_ptr = M.ptr<float>(j);
        // D holds the index of the images that contributed the maximum gradient
        uchar* D_ptr = D.ptr<uchar>(j);

        for (int k = 0; k < M.cols; k++) {
            // use the relation that |a| + |b| + |c| <= |a + b + c|
            double dx = fabs(dx0_ptr[k]) + fabs(dx1_ptr[k]) + fabs(dx2_ptr[k]);
            double dy = fabs(dy0_ptr[k]) + fabs(dy1_ptr[k]) + fabs(dy2_ptr[k]);
            // color gradient magnitude
            double grad_mag = sqrt(dx*dx + dy*dy);

            if( grad_mag > M_ptr[k]) {
                M_ptr[k] = grad_mag;
                D_ptr[k] = depth;
            }
        }
    }
}

void FocusStackWorker::coarseDepthMap(cv::Mat &dst) {

    // Use the aligned cache images. Need to read the first one to
    // get the image size
    cv::Mat src = cv::imread(cacheFilenames.at(0).toStdString());

    // Depth map output
    dst = cv::Mat::zeros(src.size(), CV_8UC1);

    // Single channel Maximum Color Gradient Magnitude
    cv::Mat G = cv::Mat::zeros(src.size(), CV_32F);
    // Single channel divergence
    cv::Mat L = cv::Mat::zeros(src.size(), CV_32F);

    // Derivative buffers
    cv::Mat dx[3], dy[3];

    // For matrix range
    double maxG, minL;

    dx[0].create(src.size(), CV_32F);
    dx[1].create(src.size(), CV_32F);
    dx[2].create(src.size(), CV_32F);

    dy[0].create(src.size(), CV_32F);
    dy[1].create(src.size(), CV_32F);
    dy[2].create(src.size(), CV_32F);

    // Gradient Magnitude
    colorDerivative(1, src, dx, dy, G, dst);
    for (int i = 1; i < cacheFilenames.count(); i++) {
        cv::Mat src = cv::imread(cacheFilenames.at(i).toStdString());
        colorDerivative(i+1, src, dx, dy, G, dst);
    }
    cv::imwrite(prefix + "_G.png", scalar2gray(G));

    // Gradient of G for edge detection
    cv::GaussianBlur(G, L, cv::Size(7,7), -1);  L.copyTo(G);
    cv::Scharr(G, dx[0], CV_32F, 1, 0);
    cv::Scharr(G, dy[0], CV_32F, 0, 1);

    cv::imwrite(prefix + "_E.png", vec2f2map(dx[0], dy[0]));

    // Second derivative for curvature filter
    cv::Laplacian(G, L, CV_32F);

    // minimum gradient and curvature filter
    cv::minMaxLoc(G, nullptr, &maxG, nullptr, nullptr);
    cv::minMaxLoc(L, &minL, nullptr, nullptr, nullptr);
    double tolG = maxG*0.1;
    double tolL = minL*0.1;

    cv::imwrite(prefix + "_D_no_filter.png", depth2map(dst));

    // Gradient Filter
    for (int j = 1; j < src.rows-1; j++) {
        float* dx_ptr = dx[0].ptr<float>(j);
        float* dy0_ptr = dy[0].ptr<float>(j-1);
        float* dy1_ptr = dy[0].ptr<float>(j);
        float* G_ptr = G.ptr<float>(j);
        float* L_ptr = L.ptr<float>(j);

        uchar* dst_ptr = dst.ptr<uchar>(j);

        for (int k = 1; k < src.cols-1; k++) {
            // backward or forward difference is required for single pixel
            // zero crossing detection. central difference always detect two
            // pixels: before and after the zero crossing.
            bool x = dx_ptr[k-1] > 0 && dx_ptr[k] < 0;
            bool y = dy0_ptr[k] > 0 && dy1_ptr[k] < 0;

            // Zero crossing in either x or y, minimin gradient and minimum
            // negative curvature
            bool rule = (x || y) && G_ptr[k] > tolG && L_ptr[k] < tolL;

            if (rule == false) {
                dst_ptr[k] = 0;
            }
        }
    }

    cv::imwrite(prefix + "_D.png", depth2map(dst));
}

void FocusStackWorker::downSampleCoarseDepthMap(cv::Mat &src, cv::Mat &dst) {
    printf("> FocusStacker::computeCoarseDenseDepthMap() ... \n"); fflush(stdout);

    int iw = src.cols;
    int ih = src.rows;

    // always downsample to a height of 300px
    int rows = 300;
    double sf = ih/rows;
    int cols = iw/sf;

    // Down sampled depth map
    dst = cv::Mat::zeros(rows, cols, CV_8UC1);

    // Sampling radius
    int r = ceil(sf/2);  if (r == 0) { r = 1; }
    // sampling window size
    int wsz = 2*r+1;

    // number of depth bins
    int dsz = cacheFilenames.count()+1;
    // depth hystogram to peek the mode (the most repeating element) in the window
    int *hyst = new int[dsz];

    printf("    Down sampling ..."); fflush(stdout);
    for (int i = r; i < rows-r; i++)
    {
        uchar *dst_ptr = dst.ptr<uchar>(i);

        for (int j = r; j < cols-r; j++)
        {
            // start row and column of the window
            int si = i * ih / rows;
            int sj = j * iw / cols;
            // discard windows that goes behond the image borders
            if (si+r >= ih || sj+r >= iw) {
                printf("    (%d) (%d)\n", si,sj);
                continue;
            }
            // Initialize the histogram
            for (int d = 0; d < dsz; d++) {
                hyst[d] = 0;
            }
            // Build the depth histogram
            for (int y = 0; y < wsz; y++) {
                uchar *W_ptr = src.ptr<uchar>(si+y);
                for (int x = 0; x < wsz; x++) {
                    // consider only depth source points (pictures index >= 1)
                    if (W_ptr[sj+x] > 0) {
                        hyst[W_ptr[sj+x]]++;
                    }
                }
            }
            // Find the most repeating depth
            int maxDepth = 0;
            int maxIndex = 0;
            for (int d = 0; d < dsz; d++) {
                // a minimum of 4 repetition are required to filter noise
                if (hyst[d] > maxDepth && hyst[d] > 4) {
                    maxDepth = hyst[d];
                    maxIndex = d;
                }
            }
            // Set it as the sampling value
            dst_ptr[j] = maxIndex;
        }
    }

    cv::imwrite(prefix + "_D_ds.png", depth2map(dst));

    printf("done\n"); fflush(stdout);
}

void FocusStackWorker::fineDepthMap(cv::Mat &src, cv::Mat &dst) {

    printf("    Building matting Laplacian ...\n"); fflush(stdout);

    // depthmap size
    int cols = src.cols;
    int rows = src.rows;

    // window radius
    int r = 1;
    // total number of unknowns
    int N = cols*rows;
    // number of elements in the window
    int Wsz = (2*r+1)*(2*r+1);
    // number of elements in the sub-matrix A
    int Asz = Wsz*Wsz;
    // lagrange multiplier
    double lambda = 100.0;

    typedef Eigen::Triplet<double> T;

    // Number of known pixels
    int knowns = 0;

    // Convert the matrix to vectors
    std::vector<T> b_list;
    std::vector<T> D_list;
    for (int i = 0, k = 0; i < rows; i++) {
        uchar *src_ptr = src.ptr<uchar>(i);
        for (int j = 0; j < cols; j++,src_ptr++,k++) {
            int nz = (*src_ptr) != 0;
            if (nz) {
                b_list.push_back(T(k,0,lambda * (*src_ptr)));
                D_list.push_back(T(k,k,lambda));
                knowns++;
            }
        }
    }

    // right hand side of equation
    Eigen::SparseMatrix<double> b(N,1);
    b.setFromTriplets(b_list.begin(), b_list.end());

    // diagonal of left hand side
    Eigen::SparseMatrix<double> D(N,N);
    D.setFromTriplets(D_list.begin(), D_list.end());

    // Total number of generated elements
    int tlen = (cols*rows - knowns) * Asz;

    // Build matting Laplacian coefficients
    std::vector<T> L_list;
    L_list.reserve(tlen);

    int *win_inds = new int[Wsz];

    for(int i = r; i < rows-r; i++) {
        uchar *src_ptr = src.ptr<uchar>(i);
        for (int j = r; j < cols-r; j++, src_ptr++) {
            // don't rocess known pixels
            if ((*src_ptr) > 0) {
                continue;
            }

            // extract window indexes
            for (int m = i-r, k = 0; m <= i+r; m++) {
                for (int n = j-r; n <= j+r; n++,k++) {
                    win_inds[k] = m * cols + n;
                }
            }

            // compute each of the A matrix contributions
            for (int m = 0; m < Wsz; m++) {
                for (int n = 0; n < Wsz; n++) {
                    int row = win_inds[m];
                    int col = win_inds[n];
                    float val = (float)(row == col) - 1.0/Wsz;
                    L_list.push_back(T(row,col,val));
                }
            }
        }
    }

    delete[] win_inds;

    Eigen::SparseMatrix<double> L(N,N);
    L.setFromTriplets(L_list.begin(), L_list.end());
    L = L + D;

#if 0
    printf("    N      = %d \n", N);
    printf("    Wsz    = %d \n", Wsz);
    printf("    Asz    = %d \n", Asz);
    printf("    knowns = %d \n", knowns);
    printf("    tlen   = %d\n", tlen);
#endif

    printf("    done\n"); fflush(stdout);

    printf("    Solving optimization problem ...\n"); fflush(stdout);
    Eigen::VectorXd X;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    solver.compute(L);
    if(solver.info() != Eigen::Success) {
      // decomposition failed
      printf("    [Error] : decomposition failed\n"); fflush(stdout);
    } else {
        X = solver.solve(b);
        if(solver.info() != Eigen::Success) {
            printf("    [Error] : solving failed\n"); fflush(stdout);
        } else {
            printf("    Solution ready\n");
        }
    }
    printf("    done\n"); fflush(stdout);

    // Convert the solution x to openCV matrix
    for (int i = 0, k = 0; i < rows; i++) {
        uchar *src_ptr = src.ptr<uchar>(i);
        for (int j = 0; j < cols; j++, k++) {
            src_ptr[j] = (uchar)round(X(k));
        }
    }
    cv::bilateralFilter(src,dst,5,10,10);

    cv::imwrite(prefix + "_D_ds_dense.png", depth2map(dst));
}

void FocusStackWorker::stack(cv::Mat &depthmap, cv::Mat &dst) {

    dst = cv::Mat(depthmap.size(), CV_8UC3);

    for (int f = 0; f < cacheFilenames.count(); f++) {
        // determine contribution of each aligned image
        cv::Mat img = cv::imread(cacheFilenames.at(f).toStdString());

        for (int i = 0; i < depthmap.rows; i++)
        {
            uchar* depthmap_ptr = depthmap.ptr<uchar>(i);

            cv::Vec3b* dst_ptr = dst.ptr<cv::Vec3b>(i);
            cv::Vec3b* img_ptr = img.ptr<cv::Vec3b>(i);

            for (int j = 0; j < depthmap.cols; j++)
            {
                // depthmap valid index start at 1, so we need to adjust it
                int idx = (int) depthmap_ptr[j] - 1;
                // clip the image index
                if (idx > (int) cacheFilenames.count()) {
                    idx = cacheFilenames.size()-1;
                } else if (idx < 0){
                    idx = 0;
                }

                if (idx == f) {
                    dst_ptr[j] = img_ptr[j];
                }
            }
        }
    }
}

void FocusStackWorker::clear() {
    groupIds.clear();
    groups.clear();
}

int FocusStackWorker::count() {
    return groups.count();
}

void FocusStackWorker::addGroup(int id, const QStringList &filenames) {
    groupIds.append(id);
    groups.append(filenames);
}

void FocusStackWorker::upSampleFineDepthMap(cv::Mat &src, cv::Mat &dst, cv::Size size) {
    cv::Mat temp;

    cv::resize(src, temp, size, 0, 0, cv::INTER_LINEAR);
    cv::bilateralFilter(temp, dst, 9, 20,20);
    cv::imwrite(prefix + "_depthmap.png", depth2map(dst));
}

void FocusStackWorker::run() {
    for (int g = 0; g < groups.count(); g++) {
        currentGroup = g;

        // List of images in the group
        QStringList filenames = groups.at(g);

        // Location to place intermediary files
        QFileInfo fileinfo(filenames.at(0));
        prefix = cacheFolder.toStdString() + "/" + fileinfo.baseName().toStdString();

        // Aligned images are cached
        cacheFilenames.clear();
        for (int i = 0; i < filenames.count(); i++) {
            QFileInfo fileinfo(filenames.at(i));
            cacheFilenames.append(QString("%1/%2_aligned.png").arg(cacheFolder).arg(fileinfo.baseName()));
        }

        // Images at different focus have different zoom (and possibly possition)
        align();

        // Buffer matrices
        cv::Mat A, B, C;

        // Build the depthmap
        coarseDepthMap(A);
        downSampleCoarseDepthMap(A, B);
        fineDepthMap(B, C);
        upSampleFineDepthMap(C, A, A.size());

        // Merge the All-in-focus image
        stack(A, B);

        std::vector<int> params;
        params.push_back(cv::IMWRITE_JPEG_QUALITY);
        params.push_back(95);

        std::string outputFilename = outputFolder.toStdString() + "/" + fileinfo.baseName().toStdString() + "_focused.jpg";
        cv::imwrite(outputFilename, B, params);

        emit imageFocused(groupIds.at(g));
    }

    emit workerFinished();
}

void FocusStackWorker::setCacheFolder(const QString &folder) {
    cacheFolder = folder;
}

void FocusStackWorker::setOutputFolder(const QString &folder) {
    outputFolder = folder;
}



