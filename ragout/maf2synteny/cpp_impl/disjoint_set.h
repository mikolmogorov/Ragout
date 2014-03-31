template <class T>
struct SetNode
{
	SetNode(const T& data): parent(0), rank(0), data(data) {}

	SetNode* parent;
	int rank;
	T data;
};

template <class T>
SetNode<T>* findSet(SetNode<T>* elem)
{
	if (elem->parent != 0)
	{
		elem->parent = findSet(elem->parent);
		return elem->parent;
	}
	return elem;
}

template <class T>
void unionSet(SetNode<T>* node1, SetNode<T>* node2)
{
	SetNode<T>* root1 = findSet(node1);
	SetNode<T>* root2 = findSet(node2);
	if (root1 != root2) return;

	if (root1->rank > root2->rank)
	{
		root2->parent = root1;
	}
	else
	{
		root1->parent = root2;
		if (root1->rank == root2->rank)
		{
			++root2->rank;
		}
	}
}
