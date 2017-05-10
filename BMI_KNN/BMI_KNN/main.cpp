#include<iostream>
#include<fstream>
#include"KNN.h"

using namespace std;

int main()
{
	myKNN knn;
	ifstream ribodata, nonribodata, experiments;
	ifstream ribonames, nonribonames;
	ribodata.open("datafile/ribo-data.txt");
	nonribodata.open("datafile/nonribo-data.txt");
	experiments.open("datafile/experiments.txt");
	ribonames.open("datafile/ribo-names.txt");
	nonribonames.open("datafile/nonribo-names.txt");

	knn.intputSamples(ribodata, ribonames, 121, RIBOSOMAL);
	knn.intputSamples(nonribodata, nonribonames, 2346, NONRIBO);
	knn.getExperiments(experiments, COLUMNS);


	knn.K = 5;
	cout << "Insert p : ";
	cin >> knn.p;

	knn.crossValidation();


	//knn.printGenes();

	ribodata.close();
	nonribodata.close();
	experiments.close();
	ribonames.close();
	nonribonames.close();
}