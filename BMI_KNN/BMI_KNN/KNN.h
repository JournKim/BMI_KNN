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
Use Euclidean distance : �� ���� ������ �� (ex : �����ﰢ�� ����)
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

// � gene A�κ��� �ٸ� gene B������ �Ÿ�, ��� index�� ����Ǿ��ִ���
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

// �� gene ������ Euclidean Distance�� ����Ѵ�.
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

	// gene�� ������ ����ϴ� �׽�Ʈ�ڵ�
	void printGenes()
	{
		for (int i = 0; i < genes.size(); i++)
		{
			printf("%5d. %s\t : %d\n", i, genes[i].name, genes[i].type);
		}
	}

	// experiments�� ������ �Է¹޴� �Լ�.
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

	// gene���� ������ �Է�.
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

	// K-nearest neighbor�� �߿� Ribosomal�� ������ p �̻��̸� return true;
	bool classification(const gene& g)
	{
		vector<geneDistance> d; // �ٸ� gene����� distance���� �����ϴ� vector.
		geneDistance tmp;

		// ��� gene data�� ���� ���ο� gene���� �Ÿ��� ���Ѵ�.
		for (int i = 0; i < genes.size(); i++)
		{
			tmp.dist = calcDist(g, genes[i]);
			tmp.idx = i;
			tmp.type = genes[i].type;
			d.push_back(tmp);
		}

		// �Ÿ��� ����� ������ ����.
		sort(d.begin(), d.end());

		int riboCount = 0;
		// K���� ����� vector�� ribosomal�� ������ ����.
		for (int i = 0; i < K; i++)
		{
			if (d[i].type == RIBOSOMAL)
			{
				riboCount++;
			}
		}

		// K�� �߿��� ribosomal�� ������ p �̻��̸� return true.
		if (((double)riboCount / K) >= p)
			return true;
		// �׷��� ������ return false.
		else
			return false;
	}

	// Cross Validation�� �� �� ���.
	bool classificationForCV(const gene& g, vector<int>& v)
	{
		vector<geneDistance> d; // �ٸ� gene����� distance���� �����ϴ� vector.
		geneDistance tmp;

		// ��� gene data�� ���� ���ο� gene���� �Ÿ��� ���Ѵ�.
		for (int i = 0; i < v.size(); i++)
		{
			tmp.dist = calcDist(g, genes[v[i]]);
			tmp.idx = i;
			tmp.type = genes[v[i]].type;
			d.push_back(tmp);
		}

		// �Ÿ��� ����� ������ ����.
		sort(d.begin(), d.end());

		int riboCount = 0;
		// K���� ����� vector�� ribosomal�� ������ ����.
		for (int i = 0; i < K; i++)
		{
			if (d[i].type == RIBOSOMAL)
			{
				riboCount++;
			}
		}

		// K�� �߿��� ribosomal�� ������ p �̻��̸� return true.
		if (((double)riboCount / K) >= p)
			return true;
		// �׷��� ������ return false.
		else
			return false;
	}


	 //�Է������� �׷��� ������ validation test.
	void crossValidation()
	{
		vector<int> riboSet[6];
		vector<int> nonriboSet[6];
		
		int TP = 0; // Ribo�ε� Ribo��� ��
		int FP = 0; // nonRibo�ε� Ribo��� ��
		int TN = 0; // nonRibo�ε� nonRibo��� ��
		int FN = 0; // Ribo�ε� nonRibo��� ��.


		// 6���� �׷����� ����.
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
				// i�� ������ �׷��� gene���� testset�� �ִ´�.
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
			// 5���� �׷��� testset�� ���� �� i��° �׷쿡 ���� test�Ѵ�.

			for (int j = 0; j < riboSet[i].size(); j++)
			{
				//cout << "ribo i : " << i << endl;
				// ribo�ε� ribo��� �Ǻ�
				if (classificationForCV(genes[riboSet[i][j]], testset))
				{
					TP++;
				}
				else // ribo�ε� nonribo��� �Ǻ�
				{
					FN++;
				}
			}
			for (int j = 0; j < nonriboSet[i].size(); j++)
			{
				//cout << "nonribo i : " << i << endl;
				//nonribo�ε� ribo��� �Ǻ�
				if (classificationForCV(genes[nonriboSet[i][j]], testset))
				{
					FPList.push_back(nonriboSet[i][j] - 120);
					FP++;
				}
				else // nonribo�ε� nonribo��� �Ǻ�
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
		double sensitivity = (double)TP / (double)(TP + FN); // ribosomal�� �߿� ribosomal�� �ǰ� �Ǵ��� ����
		double specificity = (double)TN / (double)(TN + FP); // nonribo�� �߿� nonribo�� �ǰ� �Ǵ��� ����
		double accuracy = (double)(TP + TN) / total; // ��ü �߿� �ǰ� �Ǵ��� ����.
		
		cout << "K : " << K << endl;
		cout << "p : " << p << endl;
		cout << "sensitivity : " << sensitivity << endl;
		cout << "specificity : " << specificity << endl;
		cout << "accuracy : " << accuracy << endl;

		// knn.out ������Ͽ� ���.
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