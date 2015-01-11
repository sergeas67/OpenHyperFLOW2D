/*************************************************
*   Object data loader (v.1.5.2)                 *
*   Copyright (C)  1995-2015 by Serge A. Suchkov *
*   Copyright policy: LGPL V3                    *
**************************************************/

#ifndef _obj_data_hpp_
#define _obj_data_hpp_

#ifdef __ICC
#include "intel_compatibility.h"
#endif //__ICC

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
//#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#include <iostream>
#include <iomanip>
#include <fstream>

//#define  _SAFE_ACCESS_
#include "utl/uarray.hpp"
using namespace std;

void* GetFileMapW(
                 int fd,
                 ssize_t FileSize
#ifdef _MPI
                ,int     rank = 0
#endif //_MPI
                 );

void* LoadSwapFile2D(char* FileName,
                     int x,
                     int y,
                     unsigned long DataSize,
                     int*  pf,
                     int* _fd,
                     ofstream* m_stream
#ifdef _MPI
                     ,int     rank = 0
#endif //_MPI
                    );

void* LoadSwapFile1D(char* FileName,
                     int x,
                     unsigned long DataSize,
                     int*  pf,
                     int* _fd,
                     ofstream* m_stream
#ifdef _MPI
                     ,int     rank = 0
#endif //_MPI
                    );

void CloseSwapFile(char*   SwapFileName,
                   void*   SwapData,
                   ssize_t  FileSize,
                   int     fd,
                   int     isDel
#ifdef _MPI
                   ,int     rank = 0
#endif //_MPI
                   );

enum InputDataError 
        {
         ID_NOERROR,
         ID_ERR_TIMEOUT
        };

enum DataType
        {
          DT_FLOAT,
          DT_INT,
          DT_STRING,
          DT_UNKNOWN
        };


class  XY_Table
{
 public:

 unsigned int         n;
 FP*              x;
 FP*              y;

 XY_Table(unsigned int N);
 ~XY_Table();
 FP GetX(unsigned int i);
 FP GetY(unsigned int i);
 unsigned int GetNumNodes();
};

extern  FP GetVal(XY_Table*,FP);

class Data
{

DataType dt;
char*    Name;
int      StrValSize;

public:

char*    StrVal;
FP   fVal;
int      iVal;

Data(char* name, FP val);
Data(char* name, int val);
Data(char* name, char* val);

~Data();

DataType GetDataType() {return dt;}
int      ConvertDataType(DataType _dt);
char*    GetName() {return Name;}

int operator  > (Data d);
int operator  < (Data d);
int operator == (Data d);
int operator != (Data d);

};

class Table : public XY_Table
{
 
char* Name;

public:

Table(char* name ,int N);
~Table();
char*  GetName();
void   SetName(char*);
FP GetVal(FP _x );
int    operator  > (Table);
int    operator  < (Table);
};

enum DATA_SOURCE {
                  DS_FILE,
                  DS_MEM,
            	  DS_PIPE
                 };

class InputData
{
 char*           Tmp_Buff;
 UArray<Table*>* tableArray;
 UArray<Data*>*  dataArray;
 static Table*   zeroTable;       
 int             fd;
 
 ifstream*       InputDataFile;
 ostream*        MessageStream;
 
 char*           DataSource;
 char*           DataName;
 
 int             err_no;
 int             DataSize;
 
 DATA_SOURCE     ds;

 int  GetDataFromFile(long timeout
#ifdef _MPI
                     ,int     rank = 0
#endif //_MPI
                      );
 int  GetDataFromMemory(
#ifdef _MPI
                         int     rank = 0
#endif //_MPI
                       );

public:

 InputData(char* DSource,DATA_SOURCE DS,ostream* fs,int data_size=0, long timeout=1
#ifdef _MPI
           ,int     rank=0
#endif //_MPI
           ); /* default timeout=1 sec*/
 ~InputData();
 
 int GetData(long timeout
#ifdef _MPI
            ,int rank = 0
#endif //_MPI
             );

 int         SaveAllDataAsText(char*);
 
 DATA_SOURCE GetDataSource();
 int         GetDataError();
 char*       GetDataName();
 
 ostream*    GetMessageStream();
 void        SetMessageStream(ostream*);

 int         GetIntVal(char* Name
#ifdef _MPI
                       , int rank = 0
#endif //_MPI
                       );
 FP      GetFloatVal(char* Name
#ifdef _MPI
                       , int rank = 0
#endif //_MPI
                         );
 char*       GetStringVal(char* Name
#ifdef _MPI
                       , int rank = 0
#endif //_MPI
                          );
 DataType    GetDataType(char* Name
#ifdef _MPI
                       , int rank = 0
#endif //_MPI
                         );
 FP      GetVal(char* Name, FP par
#ifdef _MPI
                       , int rank = 0
#endif //_MPI
                    );
 Table*      GetTable(char*
#ifdef _MPI
                     , int rank = 0
#endif //_MPI
                      );
 static Table*      GetZeroTable();
 int         GetTableSize(char*
#ifdef _MPI
                       , int rank = 0
#endif //_MPI
                          );
 int         EnumData(UArray<char*>*);
 int         EnumTable(UArray<char*>*);
 void        ClearInputData();
};

#endif // _obj_data_hpp_


