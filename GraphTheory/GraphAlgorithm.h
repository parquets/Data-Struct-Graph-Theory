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
		拓扑排序	参数：G 输入的图
					 SortResult 保存拓扑排序的结果
		        返回类型 void
	*/
	int v_num = G.vertexs();			//图的节点数量
	int* in_dex = new int[v_num];	//保存每个节点的入度
	std::queue<int> v_q;			//用于存储入度为0的图节点序号

	/* 获取图中每一个序号节点的入度 */
	for (int i = 1; i <= v_num; ++i)
		in_dex[i] = G.in_degree(i);

	/* 第一轮，先将所有入度为0的点入队*/
	for (int i = 1; i <= v_num; ++i)
		if (v_num[i] == 0) v_q.push(i);

	/* 当且仅当队列不为空，即仍然存在入度为0的节点，继续循环体 */
	while (!v_q.empty())
	{
		int tmp_v = v_q.front();
		v_q.pop();
		SortResult.push_back(tmp_v); // 弹出队列头，并将其保存在排序结果数组中

		for (int i = 1; i <= v_num; ++i)
		{
			if (i != tmp_v)
			{
				if (G.has_edge(tmp_v, i))
				{
					--in_dex[i];	//入度减一
					if (in_dex[i] == 0) v_q.push(i); //如果入度为0，则入队
				}
			}
		}
	}
}



template <typename VertexType>
void Dijkstra(const Graph<VertexType>& G, std::vector<VertexType>& PathList, int start_ver, int to_ver = -1)
{
	/* Dijkstra求最短路，参数说明： G 输入的图
								start_ver 开始节点
								to_ver 求最短路的目的节点 */

	int v_num = G.vertexs();	// 图中的节点数量
	Record<VertexType>* table = new Record[v_num + 1]; //用于记录节点的距离相关信息
	VertexType D;
	table[start_ver].prev_v = 0, table[start_ver].dis = 0;

	for (;;)
	{
		// 找到未知并且离start_ver最近的节点
		int v = MinDisUnknow(table, v_num);

		// 将该节点标记为已知
		table[v].is_know = true;

		//如果to_ver已知，则跳出循环
		if (v == to_ver) break;

		for (int i = 1; i <= v_num; ++i)
		{
			//访问v的邻接节点，并更新节点信息
			if (i != v && G.has_edge(v, i))
			{
				D = table[i].dis + G.get_weight(v, i);
				if (D < table[v].dis)
				{
					//如果新的距离小于原来的距离，则更新节点信息，包括
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
	/* 普利姆算法求最小生成树，参数说明： G 输入的图
									GenTree 生成树 */

	//获取图的节点数量
	int v_num = G.vertexs();

	//初始化生成树的大小为树的节点树
	GenTree.resize(v_num);

	//建立记录的表格
	Record<VertexType>* table = new Record<VertexType>[v_num];
	VertexType D;

	//最小生成树从第一个节点开始成长
	table[1].dis = 0;
	
	for (;;)
	{
		// 寻找未知的节点中和前一个节点距离最小的节点
		int v = MinDisUnknow(table, v_num);
		table[v].is_know = true;	//将该节点置为已知

		int prev_vertex = table[v].prev; //prev_vertex 表示前v的前向节点
		GenTree[prev_vertex].push_back(v);	//将边（prev_vertex，v）加入最小生成树

		for (int i = 1; i <= v_num; i++)
		{
			if (i != v && G.has_edge(i, v))
			{
				/* D为v到i的权重， 如果D小于节点当前的权重，则更新节点 */
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
	/* Kruskal求最短路径*/
	int v_num = G.vertexs();	//图的节点数
	int e_num = G.edges();		//图的边数
	int* father = new int[v_num + 1];	//并查集，父亲数组
	//edgeNode,存储起始节点，目的节点，边权重，节点列表表示图的所有边的集合
	edgeNode<VertexType>* edgeList = new edgeNode<VertexType>[e_num + 1];
	//优先队列，用于取出权重最小的边
	std::priority_queue<edgeNode<VertexType>> pqEdge;

	// 初始化生成树
	GenTree.resize(v_num);
	// 初始化并查集
	union_set_init(father, v_num + 1);
	// 初始化图的边列表
	setEdgeList(G, edgeList);

	// 将所有的边加入到优先队列中
	for (int i = 1; i <= e_num; i++)
		pqEdge.push(edgeList[i]);
	
	edgeNode<VertexType> tmp;
	while (!pqEdge.empty())
	{
		//找到边权最小的边
		tmp = pqEdge.top();
		pqEdge.pop();

		//找到最小边权对应的起点与终点的父亲
		int f_start = find(father, tmp.from_ver);
		int f_to = find(father, tmp.to_ver);

		//如果他们的父亲不同，即属于不同的集合，则将该边加入最小生成树
		if (f_start != f_to)
		{
			GenTree[f_start].push_back(f_to);
			//将两个集合合并
			Union(f_start, f_to);
		}
	}
}


template <typename VertexType>
void GraphSquare(const Graph<VertexType>& G)
{
	// 获取图的节点数
	int v_num = G.vertexs();

	// 对每个节点，都往下遍历两层
	for (int i = 1; i <= v_num; i++)
	{
		// 用于存储节点变量的队列
		std::queue<int> q;
		// 对节点i往下遍历一层，并将对应节点加入队列
		for (int j = 1; j <= v_num; j++)
		{
			if (G.has_edge(i, j))
				q.push(j);
		}
		//第一层相邻节点的数量获取
		int first_lay_size = q.size();
		
		for (int j = 0; j < first_lay_size; j++)
		{
			int t = q.front();
			q.pop();
			// 对第一层节点再往下遍历一层
			for (int k = 1; k <= v_num; k++)
			{
				//将第一次的下一层节点加入队列
				if (G.has_edge(t, k) && t != k)
					q.push(k);
			}
		}

		//每一个第二层子节点出队
		while (!q.empty())
		{
			int t = q.front();
			//如果原来没有边，则将边加入
			if (!G.has_edge(i, t)) G.add_edge(i, t);
			q.pop();
		}
	}
}

/* 使用邻接矩阵表示的图G，O(n)时间确定该图是否有汇，如果有，则返回该节点的标号 */
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

/* 确定一个图是否有圈 */
template <typename VertexType>
bool hasCircle(const Graph<VertexType>& G)
{
	// 获取图的节点数
	int v_num = G.vertexs();
	//标记位，用于标注节点是否已经访问过
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