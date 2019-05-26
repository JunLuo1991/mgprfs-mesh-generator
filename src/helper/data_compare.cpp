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
 This program compares given number of column statistics.
1) If there are only two columns, it will output the percentage (%)
 of win in each column
2) If there are more than two columns, the program will compute
 the average rank and standard deviation of each of column
 statistics.

 Column name should be given at the first row.

 e.g.
    Method1   Method2
    1.0       1.2
    1.7       1.5
    ...       ... 

*/

#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>
#include <SPL/math.hpp>

int main(int argc, char** argv)
{

  int num_cols = 0;
  int comp_rule = 0;
  double diff_lower_than = 0.0;
  double diff_greater_than = 0.0;

  namespace po = boost::program_options;
  po::options_description desc{ "Options" };
  desc.add_options()
    ("help,h", "print help information only")
    ("num-columns,n", po::value<int>(&num_cols),
      "specify number of columns in the data file")
    ("comp-rule,r", po::value<int>(&comp_rule)->default_value(0),
      "0/1 to specify the comparison rule.\n"
      "0: the biggest value ranks the highest\n"
      "1: the smallest value ranks the highest\n")
    ("diff-lower,l", po::value<double>(&diff_lower_than)->default_value(0.05),
     "which will print the number of cases that abs(col1-col2) <= l\n"
      "only works for two column\n")
    ("diff-greater,g", po::value<double>(&diff_greater_than)->default_value(0.10),
     "which will print the number of cases that abs(col1-col2) >= g\n"
      "only works for two column\n")
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

  if (vm.count("num-colums")) {
     if (num_cols <= 0) {
       std::cerr << "invalid number of columns: " << num_cols << "\n";
       return -1;
     }
  }

  if (vm.count("comp-rule")) {
    if (comp_rule != 0 && comp_rule != 1) {
      std::cerr << "invalid comparison rule: " << comp_rule << "\n";
      return -1;
    }
  }

  std::cerr << "number of columns: " << num_cols << "\n";
  std::vector<std::vector<double> > col_stats;
  std::vector<std::string> col_names;
  // Read column name from the first line
  for (int i = 0; i < num_cols; ++i) {
    std::string name;
    if (std::cin >> name) {
       col_names.push_back(name);
    } else {
      std::cerr << "read column name failed\n";
    }
  }


  // Read and store each column statistics in a vector
  col_stats.resize(num_cols);
  int col_count = 0;
  double data = 0;
  while (std::cin >> data) {
    col_stats[col_count].push_back(data);
    ++col_count;
    if (col_count >= num_cols) {
       col_count = 0;
    }
  }

  assert(col_stats.size() > 0);
  assert(col_stats[0].size() > 0);

  const int num_rows = col_stats[0].size();


  // if num_cols > 2, compare average ranking
  if (num_cols > 2) {

    std::vector<std::vector<double>> col_ranks(num_cols);

    // resize each vector for efficiency
    for (int i = 0; i < num_cols; ++i) {
      col_ranks[i].resize(num_rows);
    }

    // for each row, rank the statistics
    for (int i = 0; i < num_rows; ++i) {

      // a vector to store current row statistics
      std::vector<double> row_stats;
      for (int j = 0; j < num_cols; ++j) {
        row_stats.push_back(col_stats[j][i]);
      }

      // sort the statistics in the current row
      // so that their positions correspond to the rank
      if (comp_rule) {
        std::sort(row_stats.begin(), row_stats.end());
      } else {
        std::sort(row_stats.begin(), row_stats.end(), [](double a, double b){return a > b;});
      }

      // map the statistics with ranks
      std::map<double, int> rank_map;
      for (int j = 0; j < num_cols; ++j) {
        rank_map.insert(std::make_pair(row_stats[j], j+1));
      }

      // set the colum ranks based on the map
      for (int j = 0; j < num_cols; ++j) {
         if (rank_map.find(col_stats[j][i]) != rank_map.end()) {
           double rank = static_cast<double>(rank_map[col_stats[j][i]]);
           col_ranks[j][i] = rank;
         } else {
           std::cerr << "can not find element in row rank map\n";
           return -1;
         }
      }
    }

    // A vector to record the mean rank of each column
    std::vector<double> mean_ranks(num_cols);
    // A vector to record the standard deviation of each column
    std::vector<double> std_dev_ranks(num_cols);
    // For each column...
    for (int i = 0; i < num_cols; ++i) {
      double sum = std::accumulate(col_ranks[i].begin(), col_ranks[i].end(), 0.0);
      // Compute the mean rank
      double mean = sum / num_rows;
      double sqr_dev_sum = 0;
  
      for (auto&& x : col_ranks[i]) {
         sqr_dev_sum += SPL::sqr(x - mean);
      }

      // Compute the standard deviation
      std_dev_ranks[i] = std::sqrt(sqr_dev_sum / num_rows);
      mean_ranks[i] = mean;

    }

    std::cout.precision(2);
    // Output the name of each column
    for(auto&& name: col_names) {
      std::cout << name << " ";
    }
    std::cout << '\n';

    std::cout.precision(2);
    std::cout << std::fixed;
    // Output the mean rank of each column
    for (auto&& mean: mean_ranks) {
      std::cout << mean << " ";
    }
    std::cout << '\n';

    // Output the standard deviation of each column
    for (auto&& std_dev: std_dev_ranks) {
      std::cout << std_dev << " ";
    }
    std::cout << "\n\n";

  } else if (num_cols == 2) {   // If num_cols == 2, compare win ratio
      int count_1 = 0;
      int count_2 = 0;
      int count_draw = 0;

      std::vector<double> diff;

      // For each row...
      for (int i = 0; i < num_rows; ++i) {
        // compute the differences between the two columns
	double dif = col_stats[0][i] - col_stats[1][i];
	diff.push_back(dif);

        // count the number of rows that the first column win,
        // the second column wins, or draw
        if (comp_rule) {
          if (col_stats[0][i] < col_stats[1][i]) {
            ++count_1;
          } else if (col_stats[0][i] > col_stats[1][i]) {
            ++count_2;
          } else {
            ++count_draw;
          }
        } else {
          if (col_stats[0][i] > col_stats[1][i]) {
            ++count_1;
          } else if (col_stats[0][i] < col_stats[1][i]) {
            ++count_2;
          } else {
            ++count_draw;
          }
        }
      }

      // sort the differences of all rows
      std::sort(diff.begin(), diff.end());
      // get the meam difference
      double mean_diff = std::accumulate(diff.begin(), diff.end(), 0.0) / (diff.size());
      // get the minimum difference
      double min_diff = diff.front();
      // get the maximum difference
      double max_diff = diff.back();
      int mid_index = diff.size() / 2;
      // get the medium difference
      double medium_diff = diff[mid_index];
      double count_greater_than = 0;
      double count_lower_than = 0;
     
      // compute the minimum positive win margin of the first column
      double min_win_margin = max_diff;
      for(auto&& x : diff) {
         if(x > 0.01) {
		min_win_margin = x;
		break;
	 }
      }

      // count the number of cases the difference is lower than
      // a value or greater than a value
      for(auto&& x : diff) {
         if(SPL::absVal(x) <= diff_lower_than) ++count_lower_than;
         if(SPL::absVal(x) >= diff_greater_than) ++count_greater_than;
      }


      double win_1 = static_cast<double>(count_1 + count_draw) / static_cast<double>(num_rows);
      double win_2 = static_cast<double>(count_2 + count_draw) / static_cast<double>(num_rows);
      std::cout << "count_1 = " << count_1 << ", count_2 = " << count_2
                << ", count_draw = " << count_draw << ", num_rows = " << num_rows << "\n"
		<< "cases diff lower than " << diff_lower_than << " is: " << count_lower_than << "\n"
		<< "cases diff greater than " << diff_greater_than << " is: " << count_greater_than << "\n";
      std::cout.precision(2);
      std::cout << std::fixed;
      std::cout << "min diff: " << min_diff << "\n"
		<< "medium diff: " << medium_diff << "\n"
                << "max diff: " << max_diff << "\n"
		<< "win margins: " << min_win_margin << " ~ " << max_diff << "\n";

      for(auto&& name: col_names) {
        std::cout << name << "\t";
      }
      std::cout << "\n";
      std::cout.precision(3);
      std::cout << std::fixed;
      std::cout << win_1 << "\t" << win_2 << "\n";

  }
  return 0;
}
