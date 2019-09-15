#pragma once

#include <stdlib.h>
#include <stdio.h>

#define MAX_VERTEX_NUM		(1025)		 //限制一个图最大的节点数量
#define DEFAULT_VERTEX_NUM	(17)		 //一个默认图的大小
#define INFINTE				(2147483629) //一个大质数，表示无穷大，两节点之间没有连接
#define DEFAULT_WEIGHT		(1)			 //默认权值，为1可以表示默认无权。


/* 一个有向图带权图的邻接矩阵表示的图类 
	邻接矩阵adj_mat[i][j]表示第i个节点到第j个顶点的权值、
	节点的下标从1开始
*/
template <typename VertexType>
class Graph
{
public:
	Graph(int ver_num = DEFAULT_VERTEX_NUM);
	~Graph();

	int	 edges()   const;	//返回图边的数
	int  vertexs() const;	//返回图的节点数

	void add_edge(int from, int to, VertexType weight = DEFAULT_WEIGHT);	//向图添加边
	void delete_edge(int from, int to);	//从图中删除边

	bool has_edge(int from, int to)   const;	//判断图中from到to是否存在边
	VertexType	 get_weight(int from, int to) const;	//获取from到to的边权重
	VertexType*  get_adjust(VertexType v);	//返回节点v的临接节点

	int  in_degree (int vertex) const;	//返回节点vertex的入度
	int  out_degree(int vertex) const;	//返回节点vertex的出度

private:
	bool is_in(int ver_index);
	VertexType**	 adj_mat;		//图的邻接矩阵存储
	int  edge_num;		//图的边数量
	int  vertex_num;		//图的节点数量
};

template <typename VertexType>
bool Graph<VertexType>::is_in(int ver_index)
{
	/*判断节点ver_index是否在图中 */
	if (ver_index >= 1 && ver_index <= vertex_num) return true;
	else return false;
}

template <typename VertexType>
Graph<VertexType>::Graph(int ver_num)
{
	/* 图的构造函数，节点数为ver_num */
	if (ver_num > MAX_VERTEX_NUM)
	{
		printf("graph size exceed!\n");
		return;
	}

	adj_mat = new VertexType* [vertex_num + 1];
	for (int i = 0; i <= vertex_num; ++i) adj_mat[i] = new VertexType[vertex_num];
	edge_num = 0;
	vertex_num = vertex_num;
}

template <typename VertexType>
Graph<VertexType>::~Graph()
{
	/* 图的析构函数*/
	for (int i = 0; i <= vertex_num; ++i) delete[] adj_mat[i];
	delete[] adj_mat;
}

template <typename VertexType> int Graph<VertexType>::edges()	 const { return edge_num;   }
template <typename VertexType> int Graph<VertexType>::vertexs()  const { return vertex_num; }

template <typename VertexType>
void Graph<VertexType>::add_edge(int from, int to, VertexType weight)
{
	if (is_in(from) && is_in(to))
	{
		if (adj_mat[from][to] == INFINTE)
		{
			adj_mat[from][to] = weight;
			++edge_num;
		}
	}
}

template <typename VertexType>
void Graph<VertexType>::delete_edge(int from, int to)
{
	if (is_in(from) && is_in(to))
	{
		adj_mat[from][to] = INFINTE;
		--edge_num;
	}
}

template <typename VertexType>
bool Graph<VertexType>::has_edge(int from, int to) const
{
	if (is_in(from) && is_in(to)) return adj_mat[from][to] == INFINTE;
	else return false;
}

template <typename VertexType>
VertexType Graph<VertexType>::get_weight(int from, int to) const
{
	if (is_in(from) && is_in(to)) return adj_mat[from][to];
	else return INFINTE;
}

template <typename VertexType>
VertexType* Graph<VertexType>::get_adjust(VertexType v)
{
	return adj_mat[v];
}

template <typename VertexType>
int Graph<VertexType>::out_degree(int vertex) const
{
	if (!is_in(vertex)) return -1;
	int out_ = 0;
	for (int i = 1; i <= vertex_num; ++i)
		if (adj_mat[vertex][i] != INFINTE) ++out_;
	return out_;
}

template <typename VertexType>
int Graph<VertexType>::in_degree(int vertex) const
{
	if (!is_in(vertex)) return -1;
	int in_ = 0;
	for (int i = 1; i <= vertex_num; ++i)
		if (adj_mat[i][vertex] != INFINTE) ++in_;
	return in_;
}

typedef Graph<double> dGraph;
typedef Graph<int>	  iGraph;
