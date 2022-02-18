/*
 * matting.cpp
 *
 *  Created on: 2 jun. 2020
 *      Author: claud
 */

#include <QtCore/QtCore>
#include <opencv4/opencv2/opencv.hpp>
#include <eigen3/Eigen/Eigen>
#include "Matting.h"

using namespace cv;

typedef Eigen::Triplet<double> T;

void calcCovarInvMean(const double *win_I, double *win_mu, double *win_var);
void getLinearCoeff(const cv::Mat &s_alpha, const cv::Mat &s_I, cv::Mat &coeff, double epsilon = 1e-5);
void downSampleImage(const cv::Mat &I, cv::Mat &sI);
void upSampleAlphaUsingImg(const cv::Mat &d_alpha, const cv::Mat &d_I, const cv::Mat &I, cv::Mat &alpha);

void solveMultiLevelAlpha32FC3(const cv::Mat &I, const cv::Mat &consts_map, const cv::Mat &consts_vals, cv::Mat &alpha, int level) {

	int imageSize = I.rows * I.cols;
	if (imageSize > 200*200)
	{
		Mat d_I;
		Mat d_consts_map;
		Mat d_consts_vals;
		Mat d_alpha;

		downSampleImage(I, d_I);
		downSampleImage(consts_map, d_consts_map);
		downSampleImage(consts_vals, d_consts_vals);

		solveMultiLevelAlpha32FC3(d_I, d_consts_map, d_consts_vals, d_alpha, level+1);

#if 1
		cv::Mat matArray[] = {d_consts_map,d_consts_vals,d_alpha};
		cv::Mat out;
		cv::hconcat( matArray, 3, out );

		char title[64];
		sprintf(title, "Level: %d Downsample Images",level);
		cv::imshow(title, out);
#endif

		upSampleAlphaUsingImg(d_alpha, d_I, I, alpha);

#if 1
		//char title[64];
        sprintf(title, "Level: %d up sample Images",level);
        cv::imshow(title, alpha);
#endif

	}
	else
	{
		solveAlpha32FC3(I, consts_map, consts_vals, alpha);
	}
}

void solveAlpha32FC3(const cv::Mat &I, const cv::Mat &consts_map, const cv::Mat &consts_vals, cv::Mat &alpha) {
	int rows = I.rows;
	int cols = I.cols;

	printf("> Image size(%d,%d)\n",rows,cols);

	alpha.create(rows, cols, CV_32F);

	printf("> Building matting Laplacian ...\n"); fflush(stdout);

	// window radius
	const int r = 1;
	// window diameter
	const int d = 2*r+1;
	// total number of unknowns
	const int N = cols*rows;
	// number of elements in the window
	const int Wsz = d*d;
	// number of elements in the sub-matrix A
	const int Asz = Wsz*Wsz;
	// Lagrange multiplier
	double lambda = 100.0;
	// Number of known pixels
	int knowns = 0;

	printf("    Creating sparse right hand side and diagonal ... "); fflush(stdout);

	// Create sparse right hand side and diagonal
	std::vector<T> b_list;
	std::vector<T> D_list;

	for (int i = 0, k = 0; i < rows; i++)
	{
		const float *consts_map_ptr = consts_map.ptr<float>(i);
		const float *consts_vals_ptr = consts_vals.ptr<float>(i);

		for (int j = 0; j < cols; j++, k++)
		{
			if (consts_map_ptr[j] > 0.5)
			{
				knowns++;
				b_list.push_back(T(k,0,lambda * consts_vals_ptr[j]));
				D_list.push_back(T(k,k,lambda));
			}
		}
	}

	if (knowns >= rows*cols) {
	    alpha = consts_vals;
	    return;
	}

	// right hand side of equation
	Eigen::SparseMatrix<double> b(N,1);
	b.setFromTriplets(b_list.begin(), b_list.end());

	// diagonal of left hand side
	Eigen::SparseMatrix<double> D(N,N);
	D.setFromTriplets(D_list.begin(), D_list.end());

	// Total number of unknown elements
	int tlen = (cols*rows - knowns) * Asz;

	printf("done\n"); fflush(stdout);

	printf("    N      = %d \n", N);
	printf("    Wsz    = %d \n", Wsz);
	printf("    Asz    = %d \n", Asz);
	printf("    knowns = %d \n", knowns);
	printf("    tlen   = %d\n", tlen);

	printf("    Creating matrix L ... "); fflush(stdout);

	// Build matting Laplacian coefficients
	std::vector<T> L_list;
	L_list.reserve(tlen);

	// Image window's pixel index cache
	int win_inds[Wsz];

	// Image window's values cache
	double win_I_data[Wsz*3];
	double win_mu_data[3];
	double win_var_data[9];

	Mat win_I(Wsz,3,CV_64F,win_I_data);
	Mat win_mu(3,1,CV_64F,win_mu_data);
	Mat win_var(3,3,CV_64F,win_var_data);

	for(int i = r; i < rows-r; i++)
	{
		const float *consts_map_ptr = consts_map.ptr<float>(i);

		for (int j = r; j < cols-r; j++)
		{
			// don't process known pixels
			if (consts_map_ptr[j] > 0.5) {
				continue;
			}

			// extract window values
			for (int m = i-r, k = 0, p = 0; m <= i+r; m++) {
				const Vec3f *I_ptr = I.ptr<Vec3f>(m);
				for (int n = j-r; n <= j+r; n++, p++, k+=3) {
					win_inds[p] = m*cols + n;
					win_I_data[k+0] = I_ptr[n][0];
					win_I_data[k+1] = I_ptr[n][1];
					win_I_data[k+2] = I_ptr[n][2];
				}
			}

			calcCovarInvMean(win_I_data, win_mu_data, win_var_data);

			for (int k = 0; k < 9*3; k+=3) {
				win_I_data[k+0] -= win_mu_data[0];
				win_I_data[k+1] -= win_mu_data[1];
				win_I_data[k+2] -= win_mu_data[2];
			}

			Mat tvals = win_I*win_var*win_I.t();

			// compute each of the A matrix contributions
			for (int m = 0; m < Wsz; m++) {
				for (int n = 0; n < Wsz; n++) {
					int row = win_inds[n];
					int col = win_inds[m];
					double val = (double)(row == col) - (1.0 + tvals.at<double>(m,n)) / Wsz;
					L_list.push_back(T(row,col,val));
				}
			}
		}
	}

	Eigen::SparseMatrix<double> L(N,N);
	L.setFromTriplets(L_list.begin(), L_list.end());
	L = L + D;

	printf("done\n"); fflush(stdout);

	printf("    Solving optimization problem ... "); fflush(stdout);
	Eigen::VectorXd X;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(L);

	if(solver.info() != Eigen::Success) {
	  printf("[Error] : decomposition failed. "); fflush(stdout);
	}
	else {
		X = solver.solve(b);
		if(solver.info() != Eigen::Success){
			printf("[Error] : solving failed. "); fflush(stdout);
		}
	}
	printf("done\n"); fflush(stdout);

	printf("    Building alpha image ... "); fflush(stdout);
	// Convert the solution x to openCV matrix
	for (int i = 0, k = 0; i < rows; i++) {
		float *alpha_ptr = alpha.ptr<float>(i);
		for (int j = 0; j < cols; j++, k++) {
			alpha_ptr[j] = X(k);
		}
	}
	printf("done\n"); fflush(stdout);
}


void calcCovarInvMean(const double *win_I, double *win_mu, double *win_var) {

	double epsilon = 0.00001;

	// Compute the mean vector
	win_mu[0] = 0;
	win_mu[1] = 0;
	win_mu[2] = 0;
	for (int i = 0; i < 3*9; i += 3) {
		win_mu[0] += win_I[i+0];
		win_mu[1] += win_I[i+1];
		win_mu[2] += win_I[i+2];
	}
	win_mu[0] /= 9.0;
	win_mu[1] /= 9.0;
	win_mu[2] /= 9.0;

	// Compute the variance matrix
	double var[9] = {0,0,0,0,0,0,0,0,0};

	for(int i = 0; i < 3*9; i += 3) {
		var[0] += win_I[i+0] * win_I[i+0];
		var[1] += win_I[i+0] * win_I[i+1];
		var[2] += win_I[i+0] * win_I[i+2];
		var[4] += win_I[i+1] * win_I[i+1];
		var[5] += win_I[i+1] * win_I[i+2];
		var[8] += win_I[i+2] * win_I[i+2];
	}

	var[0] = var[0]/9.0 - win_mu[0]*win_mu[0] + epsilon/9.0;
	var[1] = var[1]/9.0 - win_mu[0]*win_mu[1];
	var[2] = var[2]/9.0 - win_mu[0]*win_mu[2];

	var[3] = var[1];
	var[4] = var[4]/9.0 - win_mu[1]*win_mu[1] + epsilon/9.0;
	var[5] = var[5]/9.0 - win_mu[1]*win_mu[2];

	var[6] = var[2];
	var[7] = var[5];
	var[8] = var[8]/9.0 - win_mu[2]*win_mu[2] + epsilon/9.0;

	// compute the inverse of the variance matrix
	double *A = var;
	double *B = win_var;

	double det = 0;
	det += A[0] * (A[4]*A[8] - A[5]*A[7]);
	det -= A[1] * (A[3]*A[8] - A[5]*A[6]);
	det += A[2] * (A[3]*A[7] - A[4]*A[6]);

	B[0] = (A[4]*A[8] - A[5]*A[7])/det;
	B[1] = (A[2]*A[7] - A[1]*A[8])/det;
	B[2] = (A[1]*A[5] - A[2]*A[4])/det;

	B[3] = B[1];
	B[4] = (A[0]*A[8] - A[2]*A[6])/det;
	B[5] = (A[2]*A[3] - A[0]*A[5])/det;

	B[6] = B[2];
	B[7] = B[5];
	B[8] = (A[0]*A[4] - A[1]*A[3])/det;
}

void getLinearCoeff(const cv::Mat &s_alpha, const cv::Mat &s_I, cv::Mat &coeff, double epsilon) {

	const int rows = s_I.rows;
	const int cols = s_I.cols;
	const int c = 3;

	// Linear coefficient 4D matrix
	coeff.create(rows, cols, CV_32FC(c+1));

	// window radius
	const int r = 1;
	// window diameter
	const int d = 2*r+1;
	// number of elements in the window
	const int Wsz = d*d;

	epsilon = sqrt(epsilon);

	// build kernel matrix G
	double G_data[(Wsz+c)*(c+1)];
	double GT_data[(Wsz+c)*(c+1)];
	double A_data[Wsz+c];

	cv::Mat G(Wsz+c, c+1, CV_64F, G_data);
	cv::Mat GT(c+1, Wsz+c, CV_64F, GT_data);
	cv::Mat A(Wsz+c, 1, CV_64F, A_data);

	for (size_t i = 0; i < sizeof(G_data)/sizeof(double); i++) {
	    G_data[i] = 0;
	}

	for (size_t k = 0; k < Wsz; k++) {
	    G_data[4*k+3] = 1.0;
	}

	G_data[4*(Wsz+0) + 0] = epsilon;
	G_data[4*(Wsz+1) + 1] = epsilon;
	G_data[4*(Wsz+2) + 2] = epsilon;

	for (size_t i = 0; i < sizeof(A_data)/sizeof(double); i++) {
	    A_data[i] = 0;
	}

	double tcoeff_data[4];
	cv::Mat tcoeff(4,1,CV_64F,tcoeff_data);

	for (int i = r; i < rows-r; i++) {
		Vec4f *coeff_ptr = coeff.ptr<Vec4f>(i);

		for (int j = r; j < cols-r; j++) {

			for (int m = i-r, k = 0, l = 0; m <= i+r; m++)
			{
				const Vec3f *s_I_ptr = s_I.ptr<Vec3f>(m);
				const float *s_alpha_ptr = s_alpha.ptr<float>(m);

				for (int n = j-r; n <= j+r; n++, k += 4, l += 1)
				{
					G_data[k+0] = s_I_ptr[n][0];
					G_data[k+1] = s_I_ptr[n][1];
					G_data[k+2] = s_I_ptr[n][2];
					A_data[l] = s_alpha_ptr[n];
				}
			}

			GT = G.t();

			tcoeff = (GT * G).inv() * GT * A;

			coeff_ptr[j][0] = tcoeff_data[0];
			coeff_ptr[j][1] = tcoeff_data[1];
			coeff_ptr[j][2] = tcoeff_data[2];
			coeff_ptr[j][3] = tcoeff_data[3];
		}
	}

	// top border
	for(int i = 0; i < r; i++)
		for(int j = 0; j < cols; j++)
			coeff.at<Vec4f>(i,j) = coeff.at<Vec4f>(r,j);
	// bottom border
	for(int i = rows-1; i > rows-1-r; i--)
		for(int j = 0; j < cols; j++)
			coeff.at<Vec4f>(i,j) = coeff.at<Vec4f>(rows-1-r,j);
	// left border
	for(int j = 0; j < r; j++)
		for(int i = 0; i < rows; i++)
			coeff.at<Vec4f>(i,j) = coeff.at<Vec4f>(i,r);
	// right border
	for(int j = cols-1; j > cols-1-r; j--)
		for(int i = 0; i < rows; i++)
			coeff.at<Vec4f>(i,j) = coeff.at<Vec4f>(i,cols-1-r);
}

void downSampleImage(const cv::Mat &I, cv::Mat &d_I) {
	resize(I, d_I, I.size()/2, 0,0, INTER_AREA );
}

void upSampleAlphaUsingImg(const cv::Mat &d_alpha, const cv::Mat &d_I, const cv::Mat &I, cv::Mat &alpha) {
	alpha.create(I.size(), CV_32F);

	Mat coeff;
	getLinearCoeff(d_alpha, d_I, coeff);

	Mat u_coeff;
	resize(coeff, u_coeff, I.size(), 0,0, INTER_LINEAR);

	const int rows = I.rows;
	const int cols = I.cols;

	for (int i = 0; i < rows; i++) {
		// guidance image
		const Vec3f *I_ptr = I.ptr<Vec3f>(i);
		// up-scaled linear matting coefficients
		Vec4f *u_coeff_ptr = u_coeff.ptr<Vec4f>(i);
		// up-scaled alpha matting
		float *alpha_ptr = alpha.ptr<float>(i);

		// coeff = (aR,aG,aB,b)
		// alpha = a*I + b
		for (int j = 0; j < cols; j++)
		{
			alpha_ptr[j]  = u_coeff_ptr[j][3];
			alpha_ptr[j] += u_coeff_ptr[j][2]*I_ptr[j][2];
			alpha_ptr[j] += u_coeff_ptr[j][1]*I_ptr[j][1];
			alpha_ptr[j] += u_coeff_ptr[j][0]*I_ptr[j][0];
		}
	}
}

