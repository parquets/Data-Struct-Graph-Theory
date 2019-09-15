#pragma once

#include <stdlib.h>
#include <stdio.h>

#define MAX_VERTEX_NUM		(1025)		 //����һ��ͼ���Ľڵ�����
#define DEFAULT_VERTEX_NUM	(17)		 //һ��Ĭ��ͼ�Ĵ�С
#define INFINTE				(2147483629) //һ������������ʾ��������ڵ�֮��û������
#define DEFAULT_WEIGHT		(1)			 //Ĭ��Ȩֵ��Ϊ1���Ա�ʾĬ����Ȩ��


/* һ������ͼ��Ȩͼ���ڽӾ����ʾ��ͼ�� 
	�ڽӾ���adj_mat[i][j]��ʾ��i���ڵ㵽��j�������Ȩֵ��
	�ڵ���±��1��ʼ
*/
template <typename VertexType>
class Graph
{
public:
	Graph(int ver_num = DEFAULT_VERTEX_NUM);
	~Graph();

	int	 edges()   const;	//����ͼ�ߵ���
	int  vertexs() const;	//����ͼ�Ľڵ���

	void add_edge(int from, int to, VertexType weight = DEFAULT_WEIGHT);	//��ͼ��ӱ�
	void delete_edge(int from, int to);	//��ͼ��ɾ����

	bool has_edge(int from, int to)   const;	//�ж�ͼ��from��to�Ƿ���ڱ�
	VertexType	 get_weight(int from, int to) const;	//��ȡfrom��to�ı�Ȩ��
	VertexType*  get_adjust(VertexType v);	//���ؽڵ�v���ٽӽڵ�

	int  in_degree (int vertex) const;	//���ؽڵ�vertex�����
	int  out_degree(int vertex) const;	//���ؽڵ�vertex�ĳ���

private:
	bool is_in(int ver_index);
	VertexType**	 adj_mat;		//ͼ���ڽӾ���洢
	int  edge_num;		//ͼ�ı�����
	int  vertex_num;		//ͼ�Ľڵ�����
};

template <typename VertexType>
bool Graph<VertexType>::is_in(int ver_index)
{
	/*�жϽڵ�ver_index�Ƿ���ͼ�� */
	if (ver_index >= 1 && ver_index <= vertex_num) return true;
	else return false;
}

template <typename VertexType>
Graph<VertexType>::Graph(int ver_num)
{
	/* ͼ�Ĺ��캯�����ڵ���Ϊver_num */
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
	/* ͼ����������*/
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
