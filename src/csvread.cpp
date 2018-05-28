/*! \file csvread.cpp
  \brief Read and store input data from file csv
*/

#include "csvread.h"

void 
csvread(std::ifstream& file,sparse_matrix& CF){

    int		col, row, tot_row;
	double	num;
    char 	virgola;
    std::string aux;

	std::getline(file,aux);	
	
	tot_row = 0;
	while(!file.eof()){
		std::getline(file,aux);
		file>>num;
		tot_row++;
	}

	CF.resize(tot_row);
	
	file.clear();
	file.seekg(0, std::ios::beg);
	
    std::getline(file,aux);
	
	col = 1;
	row = 0;

	while(row!=tot_row){

		if(col%5 == 0){
			file>>num;
			
			if(num == ','){
				std::cout<<"errore nel file di input"<<std::endl;
				break;
			}
			
			if(file.eof())
				break;
			
			CF[row][col-1] = num;
			
			col = 1;
			row++;
		}
		else{
			file>>num;
			file>>virgola;
			
			if(num == ','){
				std::cout<<"errore nel file di input"<<std::endl;
				break;
			}
			
			if(file.eof())
				break;
						
			CF[row][col-1] = num;
			col++;
		}
	}

   file.close();  
};