#pragma once
#include<vector>
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<istream>
#include<algorithm>

using std::cout;
using std::endl;
using std::cin;
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::sort;

#define COLUMNS 79

#define UNKNOWN 0
#define RIBOSOMAL 1
#define NONRIBO -1


/*
Use Euclidean distance : 각 성분 제곱의 합 (ex : 직각삼각형 빗변)
*/

class gene
{
public:
	vector<double> expLevel;
	string name;
	string description;
	int type = UNKNOWN;
	gene()
	{

	}
	gene(ifstream& in)
	{
		expLevel.resize(COLUMNS);
		for (int i = 0; i < COLUMNS; i++)
		{
			in >> expLevel[i];
		}
	}
	
};

// 어떤 gene A로부터 다른 gene B까지의 거리, 몇번 index에 저장되어있는지
class geneDistance {
public:
	double dist;
	int idx;
	int type;

	bool operator<(const geneDistance& d)
	{
		return dist < d.dist;
	}
};

// 두 gene 사이의 Euclidean Distance를 계산한다.
double calcDist(const gene& a, const gene& b)
{
	double sum = 0;
	for (int i = 0; i < COLUMNS; i++)
	{
		sum += (a.expLevel[i] - b.expLevel[i])*(a.expLevel[i] - b.expLevel[i]);
	}
	sum /= COLUMNS;
	return sum;
}

class myKNN
{
public:
	double p;
	int K;
	vector<gene> genes;
	vector<string> experiments;

	// gene의 정보를 출력하는 테스트코드
	void printGenes()
	{
		for (int i = 0; i < genes.size(); i++)
		{
			printf("%5d. %s\t : %d\n", i, genes[i].name, genes[i].type);
		}
	}

	// experiments의 정보를 입력받는 함수.
	void getExperiments(ifstream& in, int n)
	{
		char tmp[1024];
		for (int i = 0; i < n; i++)
		{
			in.getline(tmp, 500);
			experiments.push_back(tmp);
		}

		//for (int i = 0; i < experiments.size(); i++)
		//{
		//	cout<<experiments[i] << endl;
		//}
	}

	// gene들의 정보를 입력.
	void intputSamples(ifstream& inputFile, ifstream& names, int n, int type)
	{
		int number;
		char name[20];
		char desc[512];
		for (int i = 0; i < n; i++)
		{
			
			gene tmp = gene(inputFile);
			names >> number;
			names >> name;
			names.getline(desc, 500);
			tmp.name = name;
			tmp.description = desc;
			tmp.type = type;

			genes.push_back(tmp);
		}
	}

	// K-nearest neighbor들 중에 Ribosomal의 비율이 p 이상이면 return true;
	bool classification(const gene& g)
	{
		vector<geneDistance> d; // 다른 gene들과의 distance들을 저장하는 vector.
		geneDistance tmp;

		// 모든 gene data에 대해 새로운 gene과의 거리를 구한다.
		for (int i = 0; i < genes.size(); i++)
		{
			tmp.dist = calcDist(g, genes[i]);
			tmp.idx = i;
			tmp.type = genes[i].type;
			d.push_back(tmp);
		}

		// 거리가 가까운 순으로 정렬.
		sort(d.begin(), d.end());

		int riboCount = 0;
		// K개의 가까운 vector중 ribosomal의 개수를 센다.
		for (int i = 0; i < K; i++)
		{
			if (d[i].type == RIBOSOMAL)
			{
				riboCount++;
			}
		}

		// K개 중에서 ribosomal의 비율이 p 이상이면 return true.
		if (((double)riboCount / K) >= p)
			return true;
		// 그렇지 않으면 return false.
		else
			return false;
	}

	// Cross Validation을 할 때 사용.
	bool classificationForCV(const gene& g, vector<int>& v)
	{
		vector<geneDistance> d; // 다른 gene들과의 distance들을 저장하는 vector.
		geneDistance tmp;

		// 모든 gene data에 대해 새로운 gene과의 거리를 구한다.
		for (int i = 0; i < v.size(); i++)
		{
			tmp.dist = calcDist(g, genes[v[i]]);
			tmp.idx = i;
			tmp.type = genes[v[i]].type;
			d.push_back(tmp);
		}

		// 거리가 가까운 순으로 정렬.
		sort(d.begin(), d.end());

		int riboCount = 0;
		// K개의 가까운 vector중 ribosomal의 개수를 센다.
		for (int i = 0; i < K; i++)
		{
			if (d[i].type == RIBOSOMAL)
			{
				riboCount++;
			}
		}

		// K개 중에서 ribosomal의 비율이 p 이상이면 return true.
		if (((double)riboCount / K) >= p)
			return true;
		// 그렇지 않으면 return false.
		else
			return false;
	}


	 //입력파일의 그룹을 나눠서 validation test.
	void crossValidation()
	{
		vector<int> riboSet[6];
		vector<int> nonriboSet[6];
		
		int TP = 0; // Ribo인데 Ribo라고 평가
		int FP = 0; // nonRibo인데 Ribo라고 평가
		int TN = 0; // nonRibo인데 nonRibo라고 평가
		int FN = 0; // Ribo인데 nonRibo라고 평가.


		// 6개의 그룹으로 나눔.
		for (int i = 0; i < 121; i++)
		{
			riboSet[i % 6].push_back(i);
		}
		for (int i=121; i < genes.size(); i++)
		{
			nonriboSet[i % 6].push_back(i);
		}

		vector<int> testset;
		vector<int> FPList;

		for (int i = 0; i < 6; i++)
		{
			//cout << "testset " << i << endl;
			testset.clear();
			for (int j = 0; j < 6; j++)
			{
				// i를 제외한 그룹의 gene들을 testset에 넣는다.
				if (i != j)
				{
					// riboset
					for (int k = 0; k < riboSet[j].size(); k++)
					{
						testset.push_back(riboSet[j][k]);
					}
					// nonriboset
					for (int k = 0; k < nonriboSet[j].size(); k++)
					{
						testset.push_back(nonriboSet[j][k]);
					}
				}
			}

			//cout << "testset complete" << endl;
			// 5개의 그룹을 testset에 넣은 후 i번째 그룹에 대해 test한다.

			for (int j = 0; j < riboSet[i].size(); j++)
			{
				//cout << "ribo i : " << i << endl;
				// ribo인데 ribo라고 판별
				if (classificationForCV(genes[riboSet[i][j]], testset))
				{
					TP++;
				}
				else // ribo인데 nonribo라고 판별
				{
					FN++;
				}
			}
			for (int j = 0; j < nonriboSet[i].size(); j++)
			{
				//cout << "nonribo i : " << i << endl;
				//nonribo인데 ribo라고 판별
				if (classificationForCV(genes[nonriboSet[i][j]], testset))
				{
					FPList.push_back(nonriboSet[i][j] - 120);
					FP++;
				}
				else // nonribo인데 nonribo라고 판별
				{
					TN++;
				}
			}
		}
		/*cout << "TP : " << TP << endl;
		cout << "TN : " << TN << endl;
		cout << "FP : " << FP << endl;
		cout << "FN : " << FN << endl;
		*/
		double total = TP + TN + FP + FN;
		double sensitivity = (double)TP / (double)(TP + FN); // ribosomal을 중에 ribosomal로 옳게 판단한 비율
		double specificity = (double)TN / (double)(TN + FP); // nonribo들 중에 nonribo로 옳게 판단한 비율
		double accuracy = (double)(TP + TN) / total; // 전체 중에 옳게 판단한 비율.
		
		cout << "K : " << K << endl;
		cout << "p : " << p << endl;
		cout << "sensitivity : " << sensitivity << endl;
		cout << "specificity : " << specificity << endl;
		cout << "accuracy : " << accuracy << endl;

		// knn.out 출력파일에 출력.
		ofstream out;
		out.open("knn.out");
		out << "K : " << K << endl;
		out << "p : " << p << endl;
		out << "sensitivity : " << sensitivity << endl;
		out << "specificity : " << specificity << endl;
		out << "accuracy : " << accuracy << endl;

		/*for (int i = 0; i < FPList.size(); i++)
		{
			cout << FPList[i] << endl;
		}*/
	}
};