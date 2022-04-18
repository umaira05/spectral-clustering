#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"

using namespace std;

struct Pixel {
	int r = 0;
	int g = 0;
	int b = 0;

	Pixel() : r(0), g(0), b(0) {}
	Pixel(int r_in, int g_in, int b_in) : r(r_in), g(g_in), b(b_in) {}
};

double similarity(const Pixel p1, const Pixel p2);
vector<vector<Pixel>> read_image(ifstream& ifs, const int width, const int height);
void compress_image(vector<vector<Pixel>>& image, int& width, int& height, int new_width, int new_height);
Eigen::MatrixXd calculate_laplacian(const vector<vector<Pixel>>& image, const int width, const int height, bool norm);
Eigen::VectorXi k_means(Eigen::MatrixXd eigenvector_mat, const int clusters);
double sq_dist(const Eigen::VectorXd& p1, const Eigen::VectorXd& p2, const int dim);

int main(int argc, char ** argv) {
	string input_file = argv[1];
	int k = atoi(argv[2]);
	bool debug = (argc == 4 && (string) argv[3] == "--debug");

	ifstream img(input_file);
	string p3;
	int width, height, range;
	img >> p3 >> width >> height >> range;
	if (p3 != "P3" || range != 255) {
		cout << "Error with image" << endl;
		return 1;
	}

	vector<vector<Pixel>> image = read_image(img, width, height);
	img.close();
	cout << "Pixels extracted" << endl;

	if (height * width > 25 * 25) { // large images must be compressed
		compress_image(image, width, height, 25 * (width / sqrt(width * height)), 25 * (height / sqrt(width * height)));
		cout << "Image compressed";
		if (debug) {
			ofstream compressed("compressed.ppm");
			compressed << "P3\n" << width << " " << height << "\n255";
			for (int i = 0; i < width * height; i++) {
				if (i % width == 0) compressed << endl;
				Pixel* p = &image[i / width][i % width];
				compressed << p->r << " " << p->g << " " << p->b << " ";
			}
			compressed.close();
			cout << " as shown in compressed.ppm.";
		}
		cout << endl;
	}

	Eigen::MatrixXd norm_laplacian = calculate_laplacian(image, width, height, true);
	cout << "Normalized Laplacian calculated" << endl;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(norm_laplacian);
	cout << "Normalized eigenvectors calculated" << endl;
	if (debug) cout << "The eigenvalues are:" << endl << es.eigenvalues() << endl
				<< "The eigenvectors are:" << endl << es.eigenvectors() << endl;
	
	Eigen::VectorXi cluster_id = k_means(es.eigenvectors(), k);
	Pixel * palette = new Pixel[k];
	for (int i = 0; i < k; i++)
		palette[i] = Pixel(rand() % 255, rand() % 255, rand() % 255); // generate color palette

	ofstream output("output.ppm");
	output << "P3\n" << width << " " << height << "\n255";
	for (int i = 0; i < width * height; i++) {
		if (i % width == 0) output << endl;
		output << palette[cluster_id(i)].r << " "
			<< palette[cluster_id(i)].g << " "
			<< palette[cluster_id(i)].b << " ";
	}
	cout << "Split into " << k << " clusters as shown in output.ppm." << endl;
	output.close();
	delete[] palette;

	return 0;
}

double similarity(const Pixel p1, const Pixel p2) {
	double rmean = ((double) p1.r + p2.r) / 2.0;
	double dr = ((double)p2.r) - p1.r;
	double dg = ((double)p2.g) - p1.g;
	double db = ((double)p2.b) - p1.b;
	return exp(-1*sqrt((2 + rmean / 256) * dr * dr + 4 * dg * dg + (2 + (255 - rmean) / 256) * db * db));
}

vector<vector<Pixel>> read_image(ifstream& ifs, const int width, const int height) {
	vector<vector<Pixel>> image;
	for (int i = 0; i < height; i++) {
		vector<Pixel> row;
		for (int j = 0; j < width; j++) {
			int r, g, b;
			ifs >> r >> g >> b;
			row.push_back(Pixel(r, g, b));
		}
		image.push_back(row);
	}
	return image;
}

void compress_image(vector<vector<Pixel>>& image, int& width, int& height, int new_width, int new_height) {
	vector<vector<Pixel>> compressed_image;
	for (int i = 0; i < new_height; i++) {
		vector<Pixel> row;
		for (int j = 0; j < new_width; j++) {
			row.push_back(image[i * height / new_height][j * width / new_width]);
		}
		compressed_image.push_back(row);
	}
	image = compressed_image;
	height = new_height;
	width = new_width;
}

Eigen::MatrixXd calculate_laplacian(const vector<vector<Pixel>>& image, const int width, const int height, bool norm) {
	Eigen::MatrixXd laplacian = Eigen::MatrixXd::Zero(height * width, height * width); // degree - affinity
	Eigen::MatrixXd deg_recip_sqrt = Eigen::MatrixXd::Zero(height * width, height * width); // degree^(-1/2)

	for (int row = 0; row < height; row++) {
		for (int col = 0; col < width; col++) {
			if (row > 0) { // pixel above
				laplacian(row * width + col, (row - 1) * width + col) -= similarity(image[row][col], image[row - 1][col]);
				laplacian(row * width + col, row * width + col) += similarity(image[row][col], image[row - 1][col]);
				if (col > 0) { // pixel up and to the left
					laplacian(row * width + col, (row - 1) * width + col - 1) -= similarity(image[row][col], image[row - 1][col - 1]);
					laplacian(row * width + col, row * width + col) += similarity(image[row][col], image[row - 1][col - 1]);
				}
				if (col < width - 1) { // pixel up and to the right
					laplacian(row * width + col, (row - 1) * width + col + 1) -= similarity(image[row][col], image[row - 1][col + 1]);
					laplacian(row * width + col, row * width + col) += similarity(image[row][col], image[row - 1][col + 1]);
				}
			}
			if (row < height - 1) { // pixel below
				laplacian(row * width + col, (row + 1) * width + col) -= similarity(image[row][col], image[row + 1][col]);
				laplacian(row * width + col, row * width + col) += similarity(image[row][col], image[row + 1][col]);
				if (col > 0) { // pixel down and to the left
					laplacian(row * width + col, (row + 1) * width + col - 1) -= similarity(image[row][col], image[row + 1][col - 1]);
					laplacian(row * width + col, row * width + col) += similarity(image[row][col], image[row + 1][col - 1]);
				}
				if (col < width - 1) { // pixel down and to the right
					laplacian(row * width + col, (row + 1) * width + col + 1) -= similarity(image[row][col], image[row + 1][col + 1]);
					laplacian(row * width + col, row * width + col) += similarity(image[row][col], image[row + 1][col + 1]);
				}
			}
			if (col > 0) { // pixel to the left
				laplacian(row * width + col, row * width + col - 1) -= similarity(image[row][col], image[row][col - 1]);
				laplacian(row * width + col, row * width + col) += similarity(image[row][col], image[row][col - 1]);
			}
			if (col < width - 1) { // pixel to the right
				laplacian(row * width + col, row * width + col + 1) -= similarity(image[row][col], image[row][col + 1]);
				laplacian(row * width + col, row * width + col) += similarity(image[row][col], image[row][col + 1]);
			}
			if (norm)
				deg_recip_sqrt(row * width + col, row * width + col) = 1 / sqrt(laplacian(row * width + col, row * width + col));
		}
	}
	if (norm) laplacian = deg_recip_sqrt * laplacian * deg_recip_sqrt;
	return laplacian;
}

Eigen::VectorXi k_means(Eigen::MatrixXd eigenvector_mat, const int clusters) {
	for (int i = 0; i < eigenvector_mat.rows(); i++) { 
		double length = sqrt(sq_dist(Eigen::VectorXd::Zero(eigenvector_mat.cols()), eigenvector_mat.row(i), clusters));
		eigenvector_mat.row(i) = eigenvector_mat.row(i) * (1 / length);
	}

	Eigen::MatrixXd centroids(clusters, eigenvector_mat.cols());
	vector<int> centroid_nodes;
	for (int i = 0; i < clusters; i++) { // assign random centroids
		srand(time(0));
		int node_index = rand() % eigenvector_mat.rows();
		while (find(centroid_nodes.begin(), centroid_nodes.end(), node_index) != centroid_nodes.end())
			node_index = rand() % eigenvector_mat.rows();
		centroid_nodes.push_back(node_index);
		centroids.row(i) = eigenvector_mat.row(node_index);
	}
	Eigen::VectorXi node_assignments(eigenvector_mat.rows());
	Eigen::VectorXi prev_node_assignments(eigenvector_mat.rows());
	do {
		prev_node_assignments = node_assignments;
		for (int i = 0; i < eigenvector_mat.rows(); i++) { // assign points to closest centroid
			int closest_centroid = 0;
			double sq_dist_to_cent = sq_dist(eigenvector_mat.row(i), centroids.row(0), clusters);
			for (int j = 0; j < clusters; j++) {
				if (sq_dist(eigenvector_mat.row(i), centroids.row(j), clusters) < sq_dist_to_cent) {
					closest_centroid = j;
					sq_dist_to_cent = sq_dist(eigenvector_mat.row(i), centroids.row(j), clusters);
				}
			}
			node_assignments(i) = closest_centroid;
		}
		for (int i = 0; i < clusters; i++) { // move centroids to center of closest pixels
			Eigen::RowVectorXd new_centroid = Eigen::RowVectorXd::Zero(eigenvector_mat.cols());
			double num_points_assigned = 0;
			for (int j = 0; j < eigenvector_mat.rows(); j++) {
				if (i == node_assignments(j)) {
					num_points_assigned++;
					new_centroid = new_centroid + eigenvector_mat.row(j);
				}
			}
			new_centroid = new_centroid * (1/num_points_assigned);
			centroids.row(i) = new_centroid;
		}
	} while (!prev_node_assignments.isApprox(node_assignments));
	return node_assignments;
}

double sq_dist(const Eigen::VectorXd& p1, const Eigen::VectorXd& p2, const int dim) {
	double sum = 0;
	for (int i = 0; i < dim; i++) {
		sum += (p2(i) - p1(i)) * (p2(i) - p1(i));
	}
	return sum;
}