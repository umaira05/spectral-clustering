#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"

using namespace std;

struct Pixel {
	int r = 0;
	int g = 0;
	int b = 0;

	Pixel(int r_in, int g_in, int b_in) : r(r_in), g(g_in), b(b_in) {}
};

double similarity(const Pixel p1, const Pixel p2);
vector<vector<Pixel>> read_image(ifstream& ifs, const int width, const int height);
void compress_image(vector<vector<Pixel>>& image, int& width, int& height, int new_width, int new_height);
Eigen::MatrixXd calculate_laplacian(const vector<vector<Pixel>>& image, const int width, const int height, bool norm);
Eigen::VectorXd k_means(Eigen::MatrixXd eigenvector_mat, const int clusters);
double sq_dist(const Eigen::VectorXd& p1, const Eigen::VectorXd& p2, const int dim);

int main() {
	ifstream img("buck.ppm");
	string p3;
	int width, height, range;
	img >> p3 >> width >> height >> range;
	if (p3 != "P3" || range != 255) {
		std::cout << "Error with image" << std::endl;
		return 1;
	}

	vector<vector<Pixel>> image = read_image(img, width, height);
	img.close();

	cout << "Pixels extracted" << endl;

	if (height * width > 25 * 25) {
		compress_image(image, width, height, 25, 25);
		cout << "Image compressed" << endl;
	}

	Eigen::MatrixXd laplacian = calculate_laplacian(image, width, height, false);
	Eigen::MatrixXd norm_laplacian = calculate_laplacian(image, width, height, true);

	cout << "Laplacian calculated" << endl;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(laplacian);
	cout << "Eigenvectors calculated" << endl;

	//cout << "The eigenvalues are:" << endl << es.eigenvalues() << endl;
	//cout << "The eigenvectors are:" << endl << es.eigenvectors() << endl;
	//cout << "The eigenvector we care about is: " << endl << es.eigenvectors().col(1) << endl;
	ofstream output("output.ppm");
	output << "P3\n" << width << " " << height << "\n255";
	for (int i = 0; i < width * height; i++) {
		if (i % width == 0) output << endl;
		if (es.eigenvectors().col(1)(i) >= 0)
			output << "255 255 255 ";
		else
			output << "0 0 0 ";
	}
	cout << "Split into 2 sections as shown in output.ppm." << endl;
	output.close();

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_norm(norm_laplacian);
	cout << "Normalized eigenvectors calculated" << endl;
	//cout << "The eigenvectors are:" << endl << es_norm.eigenvectors() << endl;
	Eigen::VectorXd cluster_id = k_means(es_norm.eigenvectors(), 3);

	ofstream normoutput("normoutput.ppm");
	normoutput << "P3\n" << width << " " << height << "\n255";
	for (int i = 0; i < width * height; i++) {
		if (i % width == 0) { normoutput << endl; cout << endl; }
		if (cluster_id(i) == 0) {
			normoutput << "255 0 0 ";
			cout << "R ";
		}
		else if (cluster_id(i) == 1) {
			normoutput << "0 255 0 ";
				cout << "G ";
		}
		else if (cluster_id(i) == 2) {
			normoutput << "0 0 255 ";
			cout << "B ";
		}
		else
			normoutput << "0 0 0 ";
	}
	cout << "Split into 3 sections as shown in normoutput.ppm." << endl;
	normoutput.close();

	return 0;
}

double similarity(const Pixel p1, const Pixel p2) {
	double rmean = ((double) p1.r + p2.r) / 2.0;
	double dr = ((double)p2.r) - p1.r;
	double dg = ((double)p2.g) - p1.g;
	double db = ((double)p2.b) - p1.b;
	return 1 - sqrt((2 + rmean / 256) * dr * dr + 4 * dg * dg + (2 + (255 - rmean) / 256) * db * db)/764.83;
	//return 1 - sqrt(dr * dr + dg * dg + db * db) / (255 * sqrt(3));
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
	ofstream compressed("compressed.ppm");
	compressed << "P3\n" << width << " " << height << "\n255";
	for (int i = 0; i < width * height; i++) {
		if (i % width == 0) compressed << endl;
		Pixel* p = &compressed_image[i / width][i % width];
		compressed << p->r << " " << p->g << " " << p->b << " ";
	}
	compressed.close();
}

Eigen::MatrixXd calculate_laplacian(const vector<vector<Pixel>>& image, const int width, const int height, bool norm) {
	Eigen::MatrixXd laplacian = Eigen::MatrixXd::Zero(height * width, height * width); // degree - affinity
	Eigen::MatrixXd deg_recip_sqrt = Eigen::MatrixXd::Zero(height * width, height * width); // degree^(-1/2)

	for (int row = 0; row < height; row++) {
		for (int col = 0; col < width; col++) {
			if (row > 0) { // pixel above
				laplacian(row * width + col, (row - 1) * width + col) -= similarity(image[row][col], image[row - 1][col]);
				laplacian(row * width + col, row * width + col) += similarity(image[row][col], image[row - 1][col]);
			}
			if (row < height - 1) { // pixel below
				laplacian(row * width + col, (row + 1) * width + col) -= similarity(image[row][col], image[row + 1][col]);
				laplacian(row * width + col, row * width + col) += similarity(image[row][col], image[row + 1][col]);
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

Eigen::VectorXd k_means(Eigen::MatrixXd eigenvector_mat, const int clusters) {
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
		cout << node_index << ": " << centroids.row(i) << endl;
	}
	Eigen::VectorXd node_assignments(eigenvector_mat.rows()); 
	for (int repeat = 0; repeat < 20; repeat++) {
		for (int i = 0; i < eigenvector_mat.rows(); i++) { // assign points to closest centroid
			int closest_centroid = 0;
			double sq_dist_to_cent = sq_dist(eigenvector_mat.row(i), centroids.row(0), clusters);
			for (int j = 0; j < clusters; j++) {
				if (sq_dist(eigenvector_mat.row(i), centroids.row(j), clusters) < sq_dist_to_cent) {
					closest_centroid = j;
					sq_dist_to_cent = sq_dist(eigenvector_mat.row(i), centroids.row(j), clusters);
				}
			}
			cout << closest_centroid;
			node_assignments(i) = closest_centroid;
		}
		for (int i = 0; i < clusters; i++) {
			Eigen::RowVectorXd new_centroid = Eigen::RowVectorXd::Zero(eigenvector_mat.cols());
			double num_points_assigned = 0;
			for (int j = 0; j < eigenvector_mat.rows(); j++) {
				//if (repeat == 0)
					//cout << new_centroid << endl;
				if (i == node_assignments(j)) {
					num_points_assigned++;
					new_centroid = new_centroid + eigenvector_mat.row(j);
				}
			}
			new_centroid = new_centroid * (1/num_points_assigned);
			//cout << new_centroid - centroids.row(i) << endl;
			centroids.row(i) = new_centroid;
		}
		cout << endl;
	}
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
	return node_assignments;
}

double sq_dist(const Eigen::VectorXd& p1, const Eigen::VectorXd& p2, const int dim) {
	double sum = 0;
	for (int i = 0; i < dim; i++) {
		sum += (p2(i) - p1(i)) * (p2(i) - p1(i));
	}
	return sum;
}