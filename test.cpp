#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <tuple>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

std::tuple<std::vector<bool>, std::vector<bool>, std::vector<bool>, std::vector<int>, std::vector<int>> convert_rect(std::vector<double> data, double tolerance);
Eigen::VectorXd cross_correlation(const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y);
bool check_oscillation(const Eigen::VectorXd autocorrelation, int min_crossings, double corr_threshold);

int main(int argc, char* argv[]) {

	std::vector<double> data;
	std::ifstream infile(argv[1]);
	double num;

	while(infile >> num) {

		data.push_back(num);
	}

	infile.close();
	Eigen::Map<Eigen::VectorXd> x(&data[0], data.size());
	double mean = x.sum() / x.size();
	x.array() -= mean;
	Eigen::VectorXd test = cross_correlation(x, x);
	std::cout << check_oscillation(test, 6, 0.3) << "\n";
}

Eigen::VectorXd cross_correlation(const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y) {


	Eigen::VectorXd x_pad;
	Eigen::VectorXd y_pad;

	if(x.size() > y.size()) {

		x_pad = Eigen::VectorXd::Zero(x.size() * 2);
		x_pad.head(x.size()) = x;
		y_pad = Eigen::VectorXd::Zero(x.size() * 2);
		y_pad.head(y.size()) = y;
	}
	else {

		x_pad = Eigen::VectorXd::Zero(y.size() * 2);
		x_pad.head(x.size()) = x;
		y_pad = Eigen::VectorXd::Zero(y.size() * 2);
		y_pad.head(y.size()) = y;
	}

	Eigen::FFT<double> fft;
	Eigen::VectorXcd F_x(x_pad.size());
	Eigen::VectorXcd F_y(y_pad.size());
	fft.fwd(F_x, x_pad);
	fft.fwd(F_y, y_pad);
	Eigen::VectorXcd F_corr = F_x.conjugate().array() * F_y.array();
	Eigen::VectorXcd corr(F_corr.size());
	fft.inv(corr, F_corr);
	return corr.real().head(corr.size() / 2) / corr.real()[0];
}

bool check_oscillation(const Eigen::VectorXd autocorrelation, int min_crossings, double corr_threshold) {

	int n_crossings = 0;
	bool up = false;

	for(int i = 0; i < autocorrelation.size(); i += 1) {

		if(!up && abs(autocorrelation[i]) >= corr_threshold) {

			if(!up) {

				n_crossings += 1;
				std::cout << i << "\n";
			}

			up = true;
		}
		else if(abs(autocorrelation[i]) < corr_threshold) {

			up = false;
		}
	}

	return n_crossings >= min_crossings;
}
