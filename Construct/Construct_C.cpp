#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <math.h>
#include <iomanip> 
using namespace std;
using matrix_double = vector<vector<double>>;
using matrix_int = vector<vector<int>>;
istream &operator >> ( istream &in, matrix_int &M )
{
   M.clear();                                          // Make sure M is empty
   for ( string line; getline( in, line ); )           // Read one line at a time
   {
      stringstream ss( line );                         // Put that line in a temporary stream
      vector<int> row;
      for ( int e; ss >> e; row.push_back( e ) );   // Stream the elements on that line into a new row
      M.push_back( row );                              // Add the row to the matrix
   }
   return in;
}

istream &operator >> ( istream &in, matrix_double &M )
{
   M.clear();                                          // Make sure M is empty
   for ( string line; getline( in, line ); )           // Read one line at a time
   {
      stringstream ss( line );                         // Put that line in a temporary stream
      vector<double> row;
      for ( double e; ss >> e; row.push_back( e ) );   // Stream the elements on that line into a new row
      M.push_back( row );                              // Add the row to the matrix
   }
   return in;
}


ostream &operator << ( ostream &out, const matrix_double &M )
{
   for ( auto &row : M )
   {
      for ( auto e : row ) out << e << '\t';
      out << '\n';
   }
   return out;
}
int main(){
/********************************************************************
 *
 *                   Read the multi-block configuration file
 *
 * *****************************************************************/
	matrix_int config_m;
	ifstream conf_file("config_block.txt");
	conf_file >> config_m;
	//std::cout << config_m.size()<<std::endl;
	//std::cout << "here is okay"<<std::endl;
/********************************************************************
 *
 *                  C matrix generation
 *
 * ******************************************************************/
	int C_dim =0;
	int C_dim_total = 0;
	for(int i = 0; i < config_m.size(); i++){
	
		C_dim_total += config_m[i][1];
	}
//	std::cout << C_dim << std::endl;
	string st1 = "../Building_block/buidling_blk_";
	string filest = "/C_matrix";
	string ex = ".txt";
	string filename;
	stringstream ss;
	matrix_double C_mb;
	matrix_double C_sub;
	std::vector<double> inc_vc;
	ifstream inf;
	int k =0;
	int ele_count =0;
	for(int i = 0; i < config_m.size(); i++ ){

		k = config_m[i][0];
		ss << k;
		filename = "";
		filename = st1+ss.str()+filest+ex;
		//std::cout <<filename<< std::endl;
		inf.open(filename.c_str(), ios::in);
		inf >> C_sub;
		std::cout << C_sub.size()<<std::endl;
		ss.str("");
		inf.close();
		inf.clear();
		for(int j = 0; j < config_m[i][1]; j++){
			ele_count = 0;

			C_mb.push_back(inc_vc);
			while(ele_count <C_dim ){
			
				C_mb.back().push_back(0.0);
				ele_count++;
			}
			for(int kk = 0; kk < config_m[i][1]; kk++){
		
				C_mb.back().push_back(C_sub[j][kk]);
				ele_count++;
			}
			while(ele_count < C_dim_total){
			
				C_mb.back().push_back(0.0);
                                ele_count++;

			}
		}
		C_dim += config_m[i][1];

	}
	std::cout << C_sub[0][2]<< std::endl;
	ofstream outfile("C_Multi_block_CPU.txt");
	outfile << std::fixed << setprecision(16) << endl;
	outfile << C_mb;






}
