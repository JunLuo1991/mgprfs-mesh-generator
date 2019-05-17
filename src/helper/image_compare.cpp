// __START_OF_LICENSE__
// 
// Copyright (c) 2019 Jun Luo
// All rights reserved.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 3,
// or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with this program; see the file LICENSE.  If not,
// see <http://www.gnu.org/licenses/>.
// 
// __END_OF_LICENSE__


/*
 This program can be used to compare two images using the specified error metric
*/

#include <iostream>
#include <string>
#include <SPL/Array2.hpp>
#include <boost/program_options.hpp>
#include <opencv2/opencv.hpp>

// Reference from:
// https://docs.opencv.org/2.4/doc/tutorials/gpu/gpu-basics-similarity/gpu-basics-similarity.html



using namespace cv;

#if 0
enum class metric {
  psnr,
  ssim,
  luma_psnr
};
#endif

double get_psnr(const Mat& I1, const Mat& I2, int max_val)
{
    if (max_val != 255) {
 	std::cerr << "failed. This program only works for 8 bits image\n";
	std::exit(1);
    }

    Mat s1;
    absdiff(I1, I2, s1);       // |I1 - I2|
    s1.convertTo(s1, CV_32F);  // cannot make a square on 8 bits
    s1 = s1.mul(s1);           // |I1 - I2|^2

    Scalar s = sum(s1);         // sum elements per channel

    double sse = s.val[0] + s.val[1] + s.val[2]; // sum channels
    double max_value = static_cast<double>(max_val);

    if( sse <= 1e-10) // for small values return zero
        return 0;
    else
    {
        double mse = sse / static_cast<double>(I1.channels() * I1.total());
        double psnr = 10.0 * log10((max_value * max_value) / mse);
        return psnr;
    }
}

double get_ssim( const Mat& i1, const Mat& i2)
{
    const double C1 = 6.5025, C2 = 58.5225;
    /***************************** INITS **********************************/
    int d     = CV_32F;

    Mat I1, I2;
    i1.convertTo(I1, d);           // cannot calculate on one byte large values
    i2.convertTo(I2, d);

    Mat I2_2   = I2.mul(I2);        // I2^2
    Mat I1_2   = I1.mul(I1);        // I1^2
    Mat I1_I2  = I1.mul(I2);        // I1 * I2

    /*************************** END INITS **********************************/

    Mat mu1, mu2;   // PRELIMINARY COMPUTING
    GaussianBlur(I1, mu1, Size(11, 11), 1.5);
    GaussianBlur(I2, mu2, Size(11, 11), 1.5);

    Mat mu1_2   =   mu1.mul(mu1);
    Mat mu2_2   =   mu2.mul(mu2);
    Mat mu1_mu2 =   mu1.mul(mu2);

    Mat sigma1_2, sigma2_2, sigma12;

    GaussianBlur(I1_2, sigma1_2, Size(11, 11), 1.5);
    sigma1_2 -= mu1_2;

    GaussianBlur(I2_2, sigma2_2, Size(11, 11), 1.5);
    sigma2_2 -= mu2_2;

    GaussianBlur(I1_I2, sigma12, Size(11, 11), 1.5);
    sigma12 -= mu1_mu2;

    ///////////////////////////////// FORMULA ////////////////////////////////
    Mat t1, t2, t3;

    t1 = 2 * mu1_mu2 + C1;
    t2 = 2 * sigma12 + C2;
    t3 = t1.mul(t2);              // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))

    t1 = mu1_2 + mu2_2 + C1;
    t2 = sigma1_2 + sigma2_2 + C2;
    t1 = t1.mul(t2);               // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))

    Mat ssim_map;
    divide(t3, t1, ssim_map);      // ssim_map =  t3./t1;

    Scalar mssim = mean( ssim_map ); // mssim = average of ssim map
    return mssim[0];
}
  

int main(int argc, char** argv)
{
  std::string reference_image_file;
  std::string other_image_file;
  std::string metric;

  namespace po = boost::program_options;
  po::options_description desc{ "Options" };
  desc.add_options()
    ("help,h", "print help information only")
    ("reference-image,f", po::value<std::string>(&reference_image_file),
      "reference image file name")
    ("other-image,F", po::value<std::string>(&other_image_file),
      "other image file name.\n")
    ("metric,m", po::value<std::string>(&metric),
      "1. psnr ... peak signal to noise ratio"
      "2. ssim ... structure similarity index"
      "3. luma_psnr ... luma space psnr")
    ;
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (const po::error& e) {
    std::cerr << e.what() << '\n';
    return -1;
  }

  if (vm.count("help")) {
    std::cout << desc << '\n';
    return 1;
  }

  int max_val = 0;
  Mat img1, img2;

  if (reference_image_file.size()) {
     bool sgnd = false;
     std::vector<SPL::Array2<int>> comps;
     std::ifstream in(reference_image_file);
     if (SPL::decodePnm(in, comps, max_val, sgnd)) {
        return -1;
     }
     in.close();

     // Read reference image into OpenCV Mat object
     img1 = imread(reference_image_file, 1);
     if (!img1.data) {
        std::cerr << "Reference image is empty\n";
        return -1;
     }
  } else {
     std::cerr << "reference image is not specified.\n";
     return -1;
  }

  if (other_image_file.size()) {
     img2 = imread(other_image_file, 1);
     if (!img2.data) {
        std::cerr << "Other image is empty\n";
        return -1;
     }
  } else {
     std::cerr << "reference image is not specified.\n";
     return -1;
  }

  assert(img1.channels() == img2.channels());
  double result;

  if (metric == "psnr") {
     result = get_psnr(img1, img2, max_val);
  } else if (metric == "ssim") {
        if (img1.channels() == 1) {  // grayscale image
           result = get_ssim(img1, img2);
        } else if (img1.channels() == 3) {  // color image
           Mat img1_gray, img2_gray;
           cvtColor(img1, img1_gray, cv::COLOR_RGB2GRAY);
           cvtColor(img2, img2_gray, cv::COLOR_RGB2GRAY);
           result = get_ssim(img1_gray, img2_gray);
        } else {
           return -1;
        }
  } else if (metric == "luma_psnr") {
        if (img1.channels() == 1) {  // grayscale image
           result = get_psnr(img1, img2, max_val);
        } else if (img1.channels() == 3){  // color image
           Mat img1_gray, img2_gray;
           cvtColor(img1, img1_gray, cv::COLOR_RGB2GRAY);
           cvtColor(img2, img2_gray, cv::COLOR_RGB2GRAY);
           result = get_psnr(img1_gray, img2_gray, max_val);
        } else {
          return -1;
        }
        
  } else {
        std::cerr << "invalid error metric: " << metric << "\n";
        return -1;
  }
 
  std::cout.precision(6);
  std::cout << std::fixed << result << "\n"; 
  return 0;


}



