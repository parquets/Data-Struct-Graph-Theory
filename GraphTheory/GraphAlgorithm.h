#pragma once

#include <math.h>
#include <queue>
#include <algorithm>
#include <vector>
#include <list>
#include <stack>
#include <algorithm>
#include <set>
#include "Graph.h"

template <typename VertexType>
struct Record
{
	int  prev_v;
	VertexType	 dis;
	bool is_know;
	Record() :prev_v(-1), dis(INFINTE), is_know(false) {}
};

template <typename VertexType>
struct edgeNode
{
	int from_ver;
	int to_ver;
	VertexType weight;
	bool operator < (const edgeNode& obj) const
	{
		return weight < obj.weight;
	}
};

template <typename VertexType>
int MinDisUnknow(Record<VertexType>* TmpRecords, int RecordsNum)
{
	int min_index = -1, min_dis = INFINTE;
	for (int i = 0; i < RecordsNum; ++i)
	{
		if (!TmpRecords[i].is_know)
		{
			if (TmpRecords[i].dis < min_dis)
			{
				min_index = i;
				min_dis = TmpRecords[i].dis;
			}
		}
	}
}

template <typename VertexType>
void getPath(Record<VertexType>* R, int record_num, int target_ver, std::vector<VertexType>& PathList)
{
	std::stack<int> sk;
	sk.push(target_ver);
	int x = R[target_ver].prev;
	while (x != -1)
	{
		sk.push(x);
		x = R[x].prev;
	}
	while (!sk.empty())
	{
		int t = sk.top();
		PathList.push_back(t);
		sk.pop();
	}

}

template <typename VertexType>
void TopSort(const Graph<VertexType>& G, std::vector<int>& SortResult)
{
	/*
		��������	������G �����ͼ
					 SortResult ������������Ľ��
		        �������� void
	*/
	int v_num = G.vertexs();			//ͼ�Ľڵ�����
	int* in_dex = new int[v_num];	//����ÿ���ڵ�����
	std::queue<int> v_q;			//���ڴ洢���Ϊ0��ͼ�ڵ����

	/* ��ȡͼ��ÿһ����Žڵ����� */
	for (int i = 1; i <= v_num; ++i)
		in_dex[i] = G.in_degree(i);

	/* ��һ�֣��Ƚ��������Ϊ0�ĵ����*/
	for (int i = 1; i <= v_num; ++i)
		if (v_num[i] == 0) v_q.push(i);

	/* ���ҽ������в�Ϊ�գ�����Ȼ�������Ϊ0�Ľڵ㣬����ѭ���� */
	while (!v_q.empty())
	{
		int tmp_v = v_q.front();
		v_q.pop();
		SortResult.push_back(tmp_v); // ��������ͷ�������䱣����������������

		for (int i = 1; i <= v_num; ++i)
		{
			if (i != tmp_v)
			{
				if (G.has_edge(tmp_v, i))
				{
					--in_dex[i];	//��ȼ�һ
					if (in_dex[i] == 0) v_q.push(i); //������Ϊ0�������
				}
			}
		}
	}
}



template <typename VertexType>
void Dijkstra(const Graph<VertexType>& G, std::vector<VertexType>& PathList, int start_ver, int to_ver = -1)
{
	/* Dijkstra�����·������˵���� G �����ͼ
								start_ver ��ʼ�ڵ�
								to_ver �����·��Ŀ�Ľڵ� */

	int v_num = G.vertexs();	// ͼ�еĽڵ�����
	Record<VertexType>* table = new Record[v_num + 1]; //���ڼ�¼�ڵ�ľ��������Ϣ
	VertexType D;
	table[start_ver].prev_v = 0, table[start_ver].dis = 0;

	for (;;)
	{
		// �ҵ�δ֪������start_ver����Ľڵ�
		int v = MinDisUnknow(table, v_num);

		// ���ýڵ���Ϊ��֪
		table[v].is_know = true;

		//���to_ver��֪��������ѭ��
		if (v == to_ver) break;

		for (int i = 1; i <= v_num; ++i)
		{
			//����v���ڽӽڵ㣬�����½ڵ���Ϣ
			if (i != v && G.has_edge(v, i))
			{
				D = table[i].dis + G.get_weight(v, i);
				if (D < table[v].dis)
				{
					//����µľ���С��ԭ���ľ��룬����½ڵ���Ϣ������
					table[i].dis = D;
					table[i].prev_v = v;
				}
			}
		}
	}

	getPath(table, v_num, to_ver, PathList);

}

template <typename VertexType>
void Prim(const Graph<VertexType>& G, std::vector<std::vector<int>>& GenTree)
{
	/* ����ķ�㷨����С������������˵���� G �����ͼ
									GenTree ������ */

	//��ȡͼ�Ľڵ�����
	int v_num = G.vertexs();

	//��ʼ���������Ĵ�СΪ���Ľڵ���
	GenTree.resize(v_num);

	//������¼�ı��
	Record<VertexType>* table = new Record<VertexType>[v_num];
	VertexType D;

	//��С�������ӵ�һ���ڵ㿪ʼ�ɳ�
	table[1].dis = 0;
	
	for (;;)
	{
		// Ѱ��δ֪�Ľڵ��к�ǰһ���ڵ������С�Ľڵ�
		int v = MinDisUnknow(table, v_num);
		table[v].is_know = true;	//���ýڵ���Ϊ��֪

		int prev_vertex = table[v].prev; //prev_vertex ��ʾǰv��ǰ��ڵ�
		GenTree[prev_vertex].push_back(v);	//���ߣ�prev_vertex��v��������С������

		for (int i = 1; i <= v_num; i++)
		{
			if (i != v && G.has_edge(i, v))
			{
				/* DΪv��i��Ȩ�أ� ���DС�ڽڵ㵱ǰ��Ȩ�أ�����½ڵ� */
				D = G.get_weight(v, i);
				if (D < table[i].dis)
				{
					table[i].prev = v;
					table[i].dis = D;
				}
			}
		}
	}
}

void union_set_init(int* father, int size)
{
	for (int i = 0; i < size; ++i) father[i] = i;
}

int union_find(int* father, int x)
{
	while (father[x] != x)
	{
		x = father[x];
	}
	return x;
}

void Union(int* father, int x, int y)
{
	int fx = union_find(father, x);
	int fy = union_find(father, y);

	if (fx != fy) father[fx] = fy;
}


template <typename VertexType>
void setEdgeList(const Graph<VertexType>& G, edgeNode<VertexType>* List)
{
	int v_num = G.vertexs();
	int cnt = 0;
	edgeNode<VertexType> tmpEdge;
	for (int i = 1; i <= v_num; i++)
	{
		for (int j = 1; j <= v_num; j++)
		{
			if (G.has_edge(i, j) && i != j)
			{
				tmpEdge.from_ver = i;
				tmpEdge.to_ver = j;
				tmpEdge.weight = G.get_weight(i, j);
				List[cnt] = tmpEdge;
				cnt++;
			}
		}
	}
}

template <typename VertexType>
void Kruskal(const Graph<VertexType>& G, std::vector<std::vector<int>>& GenTree)
{
	/* Kruskal�����·��*/
	int v_num = G.vertexs();	//ͼ�Ľڵ���
	int e_num = G.edges();		//ͼ�ı���
	int* father = new int[v_num + 1];	//���鼯����������
	//edgeNode,�洢��ʼ�ڵ㣬Ŀ�Ľڵ㣬��Ȩ�أ��ڵ��б��ʾͼ�����бߵļ���
	edgeNode<VertexType>* edgeList = new edgeNode<VertexType>[e_num + 1];
	//���ȶ��У�����ȡ��Ȩ����С�ı�
	std::priority_queue<edgeNode<VertexType>> pqEdge;

	// ��ʼ��������
	GenTree.resize(v_num);
	// ��ʼ�����鼯
	union_set_init(father, v_num + 1);
	// ��ʼ��ͼ�ı��б�
	setEdgeList(G, edgeList);

	// �����еı߼��뵽���ȶ�����
	for (int i = 1; i <= e_num; i++)
		pqEdge.push(edgeList[i]);
	
	edgeNode<VertexType> tmp;
	while (!pqEdge.empty())
	{
		//�ҵ���Ȩ��С�ı�
		tmp = pqEdge.top();
		pqEdge.pop();

		//�ҵ���С��Ȩ��Ӧ��������յ�ĸ���
		int f_start = find(father, tmp.from_ver);
		int f_to = find(father, tmp.to_ver);

		//������ǵĸ��ײ�ͬ�������ڲ�ͬ�ļ��ϣ��򽫸ñ߼�����С������
		if (f_start != f_to)
		{
			GenTree[f_start].push_back(f_to);
			//���������Ϻϲ�
			Union(f_start, f_to);
		}
	}
}


template <typename VertexType>
void GraphSquare(const Graph<VertexType>& G)
{
	// ��ȡͼ�Ľڵ���
	int v_num = G.vertexs();

	// ��ÿ���ڵ㣬�����±�������
	for (int i = 1; i <= v_num; i++)
	{
		// ���ڴ洢�ڵ�����Ķ���
		std::queue<int> q;
		// �Խڵ�i���±���һ�㣬������Ӧ�ڵ�������
		for (int j = 1; j <= v_num; j++)
		{
			if (G.has_edge(i, j))
				q.push(j);
		}
		//��һ�����ڽڵ��������ȡ
		int first_lay_size = q.size();
		
		for (int j = 0; j < first_lay_size; j++)
		{
			int t = q.front();
			q.pop();
			// �Ե�һ��ڵ������±���һ��
			for (int k = 1; k <= v_num; k++)
			{
				//����һ�ε���һ��ڵ�������
				if (G.has_edge(t, k) && t != k)
					q.push(k);
			}
		}

		//ÿһ���ڶ����ӽڵ����
		while (!q.empty())
		{
			int t = q.front();
			//���ԭ��û�бߣ��򽫱߼���
			if (!G.has_edge(i, t)) G.add_edge(i, t);
			q.pop();
		}
	}
}

/* ʹ���ڽӾ����ʾ��ͼG��O(n)ʱ��ȷ����ͼ�Ƿ��л㣬����У��򷵻ظýڵ�ı�� */
template <typename VertexType>
int findMeetPoint(const Graph<VertexType>& G)
{
	int v_num = G.vertexs();
	int i = 1, j = 1;
	while (j <= v_num)
	{
		if (G.has_edge(i, j)) i++;
		else j++;
	}
	for (int i = 1; i <= v_num; i++)
	{
		if (G.has_edge(i, j) && i != j) return -1;
		if (G.has_edge(j, i) && i != j) return -1;
	}
	return i;
}

/* ȷ��һ��ͼ�Ƿ���Ȧ */
template <typename VertexType>
bool hasCircle(const Graph<VertexType>& G)
{
	// ��ȡͼ�Ľڵ���
	int v_num = G.vertexs();
	//���λ�����ڱ�ע�ڵ��Ƿ��Ѿ����ʹ�
	int* flag = new int[v_num + 1];
	memset(flag, 0, (v_num + 1) * sizeof(int));

	for (int i = 1; i <= v_num; i++)
	{
		if (flag[i] == 1) continue;
		std::queue<int> v_q;
		v_q.push(i);
		while (!v_q.empty())
		{
			int t = v_q.front();
			v_q.pop();
			for (int j = 1; j <= v_num; j++)
			{
				if (G.has_edge(t, j) && t != j)
				{
					if (flag[j] == 1) return true;
					else v_q.push(j);
				}
			}
		}
	}
	return false;
}

/*
template <typename VertexType>
bool isBinGraph(const Graph<VertexType>& G)
{
	int v_num = G.vertexs();
	int* flag = new int[v_num + 1];
	memset(flag, 0, (v_num) * sizeof(int));

	int m = 1;
	
	for (int i = 1; i <= v_num; i++)
	{
		if (flag[i] != 0)
		{
			if (flag[i] != m) return false;
			else continue;
		}
		std::queue<int> q;
		q.push(i);
		while (!q.empty())
		{
			int t = q.front();
			q.pop();
			flag[t] = m;
			m = -m;
			for (int j = 1; j <= v_num; j++)
			{
				if()
			}
		}
	}

}
*/