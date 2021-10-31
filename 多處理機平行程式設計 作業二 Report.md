多處理機平行程式設計 作業二 Report
===
F74086072 資訊112 呂文翔

---
## What have I done
### 1. Problem 1
I can't solve this problem.
### 2. Problem 2
#### 1. flow chart
![](https://i.imgur.com/xCDyplK.png)

#### 2. program functions
* `GenRandNumber` :
A random number generator.
Return the random integer number between [min, max].

* `Quicksort` : 
Used to sort the local list.

* `ComputePartner` :
Used to compute the partner of processor at different phase.
    * parameter:
        * phase, the phase of odd-even sort.
        * my_rank, the rank of processor.
    * Return:
        * the partner's rank.
```c=
int ComputePartner(int phase, int my_rank)
{
	int partner = 0;	
	if(phase % 2 == 0){				/* Even phase */
		if(my_rank % 2 != 0)			/* Odd rank */
			partner = my_rank - 1;
		else							/* Even rank */
			partner = my_rank + 1;
	}
	else{							/* Odd phase */
		if(my_rank % 2 != 0)			/* Odd rank */
			partner = my_rank + 1;
		else							/* Even rank */
			partner = my_rank - 1;
	}	
	if(partner == -1 || partner == comm_sz)
		partner = MPI_PROC_NULL;
	return partner;
};
```
* `Merge_low` :
In each phase of odd-even sort, the smaller rank processor of two partner processors keep the smaller keys and merge them to its local list.
    * parameter:
        * int *my_list, the local list of the processor.
        * int *recv_list, the local list receive from its partner.
        * int num_keys, the number of keys in local list.
```c=
void Merge_low(int *my_list, int *recv_list, int num_keys)
{
	int i = 0;
	int my_i = 0, recv_i = 0, tp_i = 0;
	int *tp_list = malloc(num_keys * sizeof(int));
	//keep the smaller keys from my_list and recv_list
	while( tp_i < num_keys){
		if( my_list[my_i] <= recv_list[recv_i]){
			tp_list[tp_i] = my_list[my_i];
			tp_i++;
			my_i++;
		}
		else{
			tp_list[tp_i] = recv_list[recv_i];
			tp_i++;
			recv_i++;
		}
	}
	//put tp_list to my_list
	for( i = 0; i < num_keys; i++)
		my_list[i] = tp_list[i];
};
```
* `Merge_high` :
In each phase of odd-even sort, the larger rank processor of two partner processors keep the larger keys and merge them to its local list.
    * parameter:
        * int *my_list, the local list of the processor.
        * int *recv_list, the local list receive from its partner.
        * int num_keys, the number of keys in local list.
```c=
void Merge_high(int *my_list, int *recv_list, int num_keys)
{
	int i = 0;
	int my_i = num_keys - 1, recv_i = num_keys - 1, tp_i = num_keys - 1;
	int *tp_list = malloc(num_keys * sizeof(int));
	//keep the smaller keys from my_list and recv_list
	while( tp_i >= 0){
		if( my_list[my_i] >= recv_list[recv_i]){
			tp_list[tp_i] = my_list[my_i];
			tp_i--;
			my_i--;
		}
		else{
			tp_list[tp_i] = recv_list[recv_i];
			tp_i--;
			recv_i--;
		}
	}
	//put tp_list to my_list
	for( i = 0; i < num_keys; i++)
		my_list[i] = tp_list[i];
};
```
#### 3. result screenshot
![](https://i.imgur.com/2BuRCDl.png)

## Analysis
### 1. Problem 1
I can't solve this problem.
### 2. Problem 2
* Observation:
    * With 4 processor, more n, the number of keys in global list, more execution time.
    * With 4 processor, n = 100000, need large time to execute.
* ![](https://i.imgur.com/TdG6BZW.png)
* Observation:
    * With the same number of processors, more n. more execution time.
    * With the same n, more processors, more execution time.
    * The input n = 100000 is too large to solve with 2 or 4 processors.
* Reason:
    * More processors means more communications between processors. Therefore, it takes more time to execution the problem with 16 processors than 2, 4, 8 processors when input is 1000 or 10000.
    * If the input n is bery large, it seems better to execute with more processors.
* ![](https://i.imgur.com/MXcKkxT.png)

## Difficulites
### 1. Problem 1
卡在scatterv的部分，根據觀察試驗的結果，推測是傳輸的數量 sendcounts的問題。原本是使用MPI create 新的MPI_Datatype來傳輸，但不斷嘗試還是失敗。
改成使用三個char array分別儲存blue, green, red的data傳輸，但還是失敗，不斷計算、嘗試後，還是找不到問題點。
跟改為少量的sendcounts(遠小於image大小)則可以成功傳輸，但沒找到傳整個image的大小。
### 2. problem 2
基本上沒有遇到困難的問題。
比起第一題，這題只需使用MPI_Gather, MPI_Send, MPI_Recv, 即可解決問題。