#include <mpi.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include "bmp.h"

using namespace std;

//定義平滑運算的次數
#define NSmooth 1000
#define ROOT 0

/*********************************************************/
/*變數宣告：                                             */
/*  bmpHeader    ： BMP檔的標頭                          */
/*  bmpInfo      ： BMP檔的資訊                          */
/*  **BMPSaveData： 儲存要被寫入的像素資料               */
/*  **BMPData    ： 暫時儲存要被寫入的像素資料           */
/*********************************************************/
BMPHEADER bmpHeader;
BMPINFO bmpInfo;
RGBTRIPLE **BMPSaveData = NULL;
RGBTRIPLE **BMPData = NULL;

int p = 0;
int rank = 0;
MPI_Comm comm = MPI_COMM_WORLD;
/*********************************************************/
/*函數宣告：                                             */
/*  readBMP    ： 讀取圖檔，並把像素資料儲存在BMPSaveData*/
/*  saveBMP    ： 寫入圖檔，並把像素資料BMPSaveData寫入  */
/*  swap       ： 交換二個指標                           */
/*  **alloc_memory： 動態分配一個Y * X矩陣               */
/*********************************************************/
int readBMP( char *fileName);        //read file
int saveBMP( char *fileName);        //save file
void swap(RGBTRIPLE *a, RGBTRIPLE *b);
RGBTRIPLE **alloc_memory( int Y, int X );        //allocate memory
void Array_2to1(int h, int w, RGBTRIPLE **darr, RGBTRIPLE *arr);

int main(int argc,char *argv[])
{
/*********************************************************/
/*變數宣告：                                             */
/*  *infileName  ： 讀取檔名                             */
/*  *outfileName ： 寫入檔名                             */
/*  startwtime   ： 記錄開始時間                         */
/*  endwtime     ： 記錄結束時間                         */
/*********************************************************/
	char *infileName = "input.bmp";
	char *outfileName = "output2.bmp";
	double startwtime = 0.0, endwtime = 0.0;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(comm, &p);
	
	//create the new mpi data type
	MPI_Datatype mpi_rgbtriple;
	int lengths[3] = {1, 1, 1};
	MPI_Aint disp[3] = {0, 1, 2};
	MPI_Datatype types[3] = {MPI_BYTE, MPI_BYTE, MPI_BYTE};

	MPI_Type_create_struct(3, lengths, disp, types, &mpi_rgbtriple);
	MPI_Type_commit(&mpi_rgbtriple);
	MPI_Comm_rank(comm, &rank);	
	
	//讀取檔案
	if( rank == ROOT){
		if ( readBMP( infileName) )
			cout << "Read file successfully!!" << endl;
		else
			cout << "Read file fails!!" << endl;
	}

	//記錄開始時間
	MPI_Barrier(comm);
	startwtime = MPI_Wtime();

	//動態分配記憶體給暫存空間
	RGBTRIPLE *BMPSaveData_local = new RGBTRIPLE[(bmpInfo.biHeight/p + 2)* bmpInfo.biWidth];	//local data, recv_buf
	RGBTRIPLE *BMPData_local = new RGBTRIPLE[(bmpInfo.biHeight/p + 2)* bmpInfo.biWidth];		//temp

	RGBTRIPLE *Data_send_1D =  new RGBTRIPLE[bmpInfo.biHeight * bmpInfo.biWidth];	//send_buf
	Array_2to1( bmpInfo.biHeight, bmpInfo.biWidth, BMPSaveData, Data_send_1D);

	//calculate sendcounts and displs
	int *sendcounts = new int[p];
	int *displs = new int[p];
	int pos = 0;
	int amount = 0;
	if( rank == 0){
		amount = bmpInfo.biHeight * bmpInfo.biWidth / p;
		for( int i = 0; i < p; i++){
			sendcounts[i] = amount;
			displs[i] = pos;
			pos += sendcounts[i];
		}
	}

	//Bcast
	MPI_Bcast(&amount, 1, MPI_INT, ROOT, comm);
	MPI_Bcast(sendcounts, p, MPI_INT, ROOT, comm);
	MPI_Bcast(displs, p, MPI_INT, ROOT, comm);
		
	//scatterd data
	MPI_Scatterv( Data_send_1D, sendcounts, displs, mpi_rgbtriple, BMPSaveData_local, amount, mpi_rgbtriple, ROOT, comm);

	cout << "rank="<< rank <<" "<<(int)Data_send_1D[100*rank].rgbBlue <<endl;


	int tp = 0;
if(tp == 1){
	//進行多次的平滑運算:
	for(int count = 0; count < NSmooth ; count ++){
		//把像素資料與暫存指標做交換
		swap(BMPSaveData,BMPData);
		//進行平滑運算
		for(int i = 0; i<bmpInfo.biHeight ; i++)
			for(int j =0; j<bmpInfo.biWidth ; j++){
				/*********************************************************/
				/*設定上下左右像素的位置                                 */
				/*********************************************************/
				int Top = i>0 ? i-1 : bmpInfo.biHeight-1;
				int Down = i<bmpInfo.biHeight-1 ? i+1 : 0;
				int Left = j>0 ? j-1 : bmpInfo.biWidth-1;
				int Right = j<bmpInfo.biWidth-1 ? j+1 : 0;
				/*********************************************************/
				/*與上下左右像素做平均，並四捨五入                       */
				/*********************************************************/
				BMPSaveData[i][j].rgbBlue =  (double) (BMPData[i][j].rgbBlue+BMPData[Top][j].rgbBlue+BMPData[Top][Left].rgbBlue+BMPData[Top][Right].rgbBlue+BMPData[Down][j].rgbBlue+BMPData[Down][Left].rgbBlue+BMPData[Down][Right].rgbBlue+BMPData[i][Left].rgbBlue+BMPData[i][Right].rgbBlue)/9+0.5;
				BMPSaveData[i][j].rgbGreen =  (double) (BMPData[i][j].rgbGreen+BMPData[Top][j].rgbGreen+BMPData[Top][Left].rgbGreen+BMPData[Top][Right].rgbGreen+BMPData[Down][j].rgbGreen+BMPData[Down][Left].rgbGreen+BMPData[Down][Right].rgbGreen+BMPData[i][Left].rgbGreen+BMPData[i][Right].rgbGreen)/9+0.5;
				BMPSaveData[i][j].rgbRed =  (double) (BMPData[i][j].rgbRed+BMPData[Top][j].rgbRed+BMPData[Top][Left].rgbRed+BMPData[Top][Right].rgbRed+BMPData[Down][j].rgbRed+BMPData[Down][Left].rgbRed+BMPData[Down][Right].rgbRed+BMPData[i][Left].rgbRed+BMPData[i][Right].rgbRed)/9+0.5;
			}
	}
}
	

	//得到結束時間，並印出執行時間
	MPI_Barrier(comm);
	endwtime = MPI_Wtime();
	if( rank == 0){ 
		cout << "The execution time = "<< endwtime-startwtime <<endl ;
	}

 	//寫入檔案
	if(rank == 0){
		if ( saveBMP( outfileName ) )
			cout << "Save file successfully!!" << endl;
		else
			cout << "Save file fails!!" << endl;
	}

	free(BMPData[0]);
 	free(BMPSaveData[0]);
 	free(BMPData);
 	free(BMPSaveData);
	//free(BMPData_local[0]);
 	//free(BMPSaveData_local[0]);
	//free(BMPData_local);
 	//free(BMPSaveData_local);
	delete sendcounts;
	delete displs;
	MPI_Type_free(&mpi_rgbtriple);
 	MPI_Finalize();

    return 0;
}

/*********************************************************/
/* 讀取圖檔                                              */
/*********************************************************/
int readBMP(char *fileName)
{
	//建立輸入檔案物件
        ifstream bmpFile( fileName, ios::in | ios::binary );

        //檔案無法開啟
        if ( !bmpFile ){
                cout << "It can't open file!!" << endl;
                return 0;
        }

        //讀取BMP圖檔的標頭資料
    	bmpFile.read( ( char* ) &bmpHeader, sizeof( BMPHEADER ) );

        //判決是否為BMP圖檔
        if( bmpHeader.bfType != 0x4d42 ){
                cout << "This file is not .BMP!!" << endl ;
                return 0;
        }

        //讀取BMP的資訊
        bmpFile.read( ( char* ) &bmpInfo, sizeof( BMPINFO ) );

        //判斷位元深度是否為24 bits
        if ( bmpInfo.biBitCount != 24 ){
                cout << "The file is not 24 bits!!" << endl;
                return 0;
        }

        //修正圖片的寬度為4的倍數
        while( bmpInfo.biWidth % 4 != 0 )
        	bmpInfo.biWidth++;

        //動態分配記憶體
        BMPSaveData = alloc_memory( bmpInfo.biHeight, bmpInfo.biWidth);

        //讀取像素資料
    	//for(int i = 0; i < bmpInfo.biHeight; i++)
        //	bmpFile.read( (char* )BMPSaveData[i], bmpInfo.biWidth*sizeof(RGBTRIPLE));
	    bmpFile.read( (char* )BMPSaveData[0], bmpInfo.biWidth*sizeof(RGBTRIPLE)*bmpInfo.biHeight);

        //關閉檔案
        bmpFile.close();

        return 1;

}
/*********************************************************/
/* 儲存圖檔                                              */
/*********************************************************/
int saveBMP( char *fileName)
{
 	//判決是否為BMP圖檔
        if( bmpHeader.bfType != 0x4d42 ){
                cout << "This file is not .BMP!!" << endl ;
                return 0;
        }

 	//建立輸出檔案物件
        ofstream newFile( fileName,  ios:: out | ios::binary );

        //檔案無法建立
        if ( !newFile ){
                cout << "The File can't create!!" << endl;
                return 0;
        }

        //寫入BMP圖檔的標頭資料
        newFile.write( ( char* )&bmpHeader, sizeof( BMPHEADER ) );

	//寫入BMP的資訊
        newFile.write( ( char* )&bmpInfo, sizeof( BMPINFO ) );

        //寫入像素資料
        //for( int i = 0; i < bmpInfo.biHeight; i++ )
        //        newFile.write( ( char* )BMPSaveData[i], bmpInfo.biWidth*sizeof(RGBTRIPLE) );
        newFile.write( ( char* )BMPSaveData[0], bmpInfo.biWidth*sizeof(RGBTRIPLE)*bmpInfo.biHeight );

        //寫入檔案
        newFile.close();

        return 1;

}


/*********************************************************/
/* 分配記憶體：回傳為Y*X的矩陣                           */
/*********************************************************/
RGBTRIPLE **alloc_memory(int Y, int X )
{
	//建立長度為Y的指標陣列
        RGBTRIPLE **temp = new RGBTRIPLE *[ Y ];
	    RGBTRIPLE *temp2 = new RGBTRIPLE [ Y * X ];
        memset( temp, 0, sizeof( RGBTRIPLE ) * Y);
        memset( temp2, 0, sizeof( RGBTRIPLE ) * Y * X );

	//對每個指標陣列裡的指標宣告一個長度為X的陣列
        for( int i = 0; i < Y; i++){
                temp[ i ] = &temp2[i*X];
        }

        return temp;

}
/*********************************************************/
/* 交換二個指標                                          */
/*********************************************************/
void swap(RGBTRIPLE *a, RGBTRIPLE *b)
{
	RGBTRIPLE *temp;
	temp = a;
	a = b;
	b = temp;
}
void Array_2to1(int h, int w, RGBTRIPLE **darr, RGBTRIPLE *arr)
{
	for(int i = 0; i < h; i++){
		for(int j = 0; j < w; j++){
			arr[i*w + j].rgbBlue = darr[i][j].rgbBlue;
			arr[i*w + j].rgbGreen = darr[i][j].rgbGreen;
			arr[i*w + j].rgbRed = darr[i][j].rgbRed;
		}
	}
}
