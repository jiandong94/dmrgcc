#include "util/genereal.h"

void QuickSort(double* array, int* index, int left, int right)
{
	if(left < right)
	{
		int i = left, j = right;
		double array_privot = array[i];
		int index_private = index[i];
		while(i < j)
		{
			while(i<j && array[j]<=array_privot)
				j--;
			if(i < j)
			{
				array[i] = array[j];
				index[i] = index[j];
				i++;
			}
			while(i<j && array[i]>=array_privot)
				i++;
			if(i < j)
			{
				array[j] = array[i];
				index[j] = index[i];
				j--;
			}
		}
		array[i] = array_privot;
		index[i] = index_private;
		QuickSort(array, index, left, i-1);
		QuickSort(array, index, i+1, right);
	}
}
