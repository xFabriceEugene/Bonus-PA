#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int extraMemoryAllocated;

void *Alloc(size_t sz)
{
	extraMemoryAllocated += sz;
	size_t* ret = malloc(sizeof(size_t) + sz);
	*ret = sz;
	printf("Extra memory allocated, size: %ld\n", sz);
	return &ret[1];
}

void DeAlloc(void* ptr)
{
	size_t* pSz = (size_t*)ptr - 1;
	extraMemoryAllocated -= *pSz;
	printf("Extra memory deallocated, size: %ld\n", *pSz);
	free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
	return ((size_t*)ptr)[-1];
}


void heapify(int arr[], int n, int i) {
    // Find largest among root, left child and right child
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
  
    if (left < n && arr[left] > arr[largest])
      largest = left;
  
    if (right < n && arr[right] > arr[largest])
      largest = right;
  
    // Swap and continue heapifying if root is not largest
    if (largest != i) {
      	int tmp = arr[i];
    	arr[i] = arr[largest];
    	arr[largest] = tmp;

      	heapify(arr, n, largest);
    }
  }

// implements heap sort
// extraMemoryAllocated counts bytes of memory allocated
void heapSort(int arr[], int N) {
    for (int i = N / 2 - 1; i >= 0; i--) {
        heapify(arr, N, i);
    }

    for (int i = N - 1; i > 0; i--) {
        int tmp = arr[0];
        arr[0] = arr[i];
        arr[i] = tmp;

        heapify(arr, i, 0);
    }
}

// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated
void mergeSort(int pData[], int l, int r){
	int m = 0;
	if(l < r){
		m = (r + l)/2;
		mergeSort(pData,l,m);
		mergeSort(pData,m + 1, r);

		int i, j, k;
		int tmpLeftSize = m - l + 1;
		int tmpRightSize = r - m;
		

		int *tmpLeftArray = Alloc(sizeof(int) * tmpLeftSize);
		int *tmpRightArray = Alloc(sizeof(int) * tmpRightSize);


		for (i = 0; i < tmpLeftSize; i++){
			tmpLeftArray[i] = pData[l + i];
		}
		for (j = 0; j < tmpRightSize; j++){
			tmpRightArray[j] = pData[m + 1+ j];
		}
		
		i = 0;
		j = 0; 
		k = l; 


		while (i < tmpLeftSize && j < tmpRightSize)
		{
			if (tmpLeftArray[i] <= tmpRightArray[j])
			{
				pData[k] = tmpLeftArray[i];
				i++;
			}
			else
			{
				pData[k] = tmpRightArray[j];
				j++;
			}
			k++;
		}


		while (i < tmpLeftSize)
		{
			pData[k] = tmpLeftArray[i];
			i++;
			k++;
		}


		while (j < tmpRightSize)
		{
			pData[k] = tmpRightArray[j];
			j++;
			k++;
		}


		DeAlloc(tmpLeftArray);
		DeAlloc(tmpRightArray);
	}

}

// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
void insertionSort(int* pData, int n)
{
	int i, key, j;
    for (i = 1; i < n; i++) {
        key = pData[i];
        j = i - 1;
 
        /* Move elements of arr[0..i-1], that are
          greater than key, to one position ahead
          of their current position */
        while (j >= 0 && pData[j] > key) {
            pData[j + 1] = pData[j];
            j = j - 1;
        }
        pData[j + 1] = key;
    }
}

// implement bubble sort
// extraMemoryAllocated counts bytes of extra memory allocated
void bubbleSort(int* pData, int n)
{
    while(n > 0){
        int tmpIndex = 0;
        for(int b = 1;b < n;++b){

            if(pData[tmpIndex]>pData[b]){
				
				
				int tmp = pData[tmpIndex];
    			pData[tmpIndex] = pData[b];
    			pData[b] = tmp;
                tmpIndex = b;
            }
            else{
                tmpIndex = b;
            }
        }
       
        --n;
	}
}

// implement selection sort
// extraMemoryAllocated counts bytes of extra memory allocated
void selectionSort(int* pData, int n)
{
    int i = n - 1;
    int baseIndex = 0;
    while(i >= 0){
        int tmpIndex = 0;
        int minIndex = baseIndex;
        for(int b = baseIndex;b < n ;++b){
            if(pData[tmpIndex]>pData[b] && pData[b]<pData[minIndex]){
                minIndex = b;
                tmpIndex = b;
            }
            else{
                tmpIndex = b;
            }
        }
        if(minIndex != baseIndex){
            int tmp = pData[minIndex];
    		pData[minIndex] = pData[baseIndex];
    		pData[baseIndex] = tmp;
            
        }
        ++baseIndex;
        --i;
    }
}

// parses input file to an integer array
int parseData(char *inputFileName, int **ppData)
{
	FILE* inFile = fopen(inputFileName,"r");
	int dataSz = 0;
	*ppData = NULL;
	int i, n, *data;
	
	if (inFile)
	{
		fscanf(inFile,"%d\n",&dataSz);
		*ppData = (int *)Alloc(sizeof(int) * dataSz);
		// Implement parse data block
		if (*ppData == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(-1);
		}
		for (i = 0; i < dataSz;++i)
		{
			fscanf(inFile, "%d ",&n);
			data = *ppData + i;
			*data = n;
		}

		fclose(inFile);
	}
	
	return dataSz;
}

// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
	int i, sz = dataSz - 100;
	printf("\tData:\n\t");
	for (i=0;i<100;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\t");
	
	for (i=sz;i<dataSz;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\n");
}

int main(void)
{
	clock_t start, end;
	int i;
    double cpu_time_used;
	char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};
	
	for (i=0;i<3;++i)
	{
		int *pDataSrc, *pDataCopy;
		int dataSz = parseData(fileNames[i], &pDataSrc);
		
		if (dataSz <= 0)
			continue;
		
		pDataCopy = (int *)Alloc(sizeof(int)*dataSz);
	
		printf("---------------------------\n");
		printf("Dataset Size : %d\n",dataSz);
		printf("---------------------------\n");
		
		printf("Selection Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		selectionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Insertion Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		insertionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Bubble Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		bubbleSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Merge Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		mergeSort(pDataCopy, 0, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
        printf("Heap Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		heapSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
	
		DeAlloc(pDataCopy);
		DeAlloc(pDataSrc);
	}
	
}
