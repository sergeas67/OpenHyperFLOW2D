/*************************************************
*   Object data loader (v.1.5.2)                 *
*   Copyright (C)  1995-2014 by Serge A. Suchkov *
*   Copyright policy: LGPL V3                    *
**************************************************/
#include "obj_data/obj_data.hpp"

#ifndef  max
    #define max(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef  min
    #define min(a,b) (((a)<(b))?(a):(b))
#endif

#ifdef _NO_MMAP_SHARED_
    #define WRITE_LARGE_FILE
#else
    #ifdef _NO_MMAP_
        #define WRITE_LARGE_FILE
    #else
        #define MSYNC_LARGE_FILE
    #endif //_NO_MMAP_
#endif // _NO_MMAP_SHARED_

static const char _STRING_TYPE[]="_STRING_";
static const char _FLOAT_TYPE[]="_FLOAT_";
static const char _INT_TYPE[]="_INT_";
static char  Message[2048];

/* -- Start file func --- */
void* GetFileMapW(int fd,
                  ssize_t FileSize
#ifdef _MPI
                  ,int rank
#endif //_MPI
                 ) {
    void*   ObjectPtr;
#ifdef _NO_MMAP_SHARED_
    ObjectPtr = (void*)mmap(0,FileSize,PROT_READ|PROT_WRITE,MAP_PRIVATE,fd,0);
#else
#ifdef WRITE_LARGE_FILE
#warning "Read file > 2GB..."
    ssize_t max_write = 1024L*1024L*1024L;
    ssize_t one_write = 0L;
    ssize_t len  = 0L;
    off_t  off = 0L;
    ObjectPtr = (void*)malloc(FileSize);
    char*  TmpPtr=(char*)ObjectPtr;
    lseek(fd,0,SEEK_SET);

    if ( FileSize > max_write ) {
        for ( off = 0L,one_write = max_write; len < FileSize; off += max_write ) {
            len += pread(fd,TmpPtr+off,one_write,off);
            if ( FileSize - len < max_write )
                one_write = FileSize - len;
        }
        if ( len != FileSize )
#ifdef _MPI
            if ( rank == 0 )
#endif //_MPI
                cout << "Error: len(" << len << ") != FileSize(" << FileSize << ") " << endl << flush;
    } else {
        len = pread(fd,ObjectPtr,FileSize,0L);
    }
#else
    ObjectPtr = (void*)mmap(0,FileSize,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
#endif // WRITE_LARGE_FILE
#endif // _NO_MMAP_SHARED_
    if ( ObjectPtr == MAP_FAILED ) ObjectPtr = NULL;
    return(ObjectPtr);
}

void CloseSwapFile(char*   SwapFileName,
                   void*   SwapData,
                   ssize_t  FileSize,
                   int     fd,
                   int     isDel
#ifdef _MPI
                   ,int rank
#endif //_MPI
                  ) {
#ifdef WRITE_LARGE_FILE
#warning "Write file > 2Gb"    
    ssize_t max_write = 1024L*1024L*1024L;
    ssize_t one_write = 0L;
    ssize_t len  = 0L;
    off_t  off = 0L;
    char*  TmpPtr=(char*)SwapData;
    if ( FileSize > max_write ) {
        for ( off = 0L,one_write = max_write; len < FileSize; off += max_write ) {
            len += pwrite64(fd,TmpPtr+off,one_write,off);
            if ( FileSize - len < max_write )
                one_write = FileSize - len;
        }
        if ( len != FileSize )
#ifdef _MPI
            if ( rank == 0 )
#endif //_MPI
                cout << "Error: len(" << len << ") != FileSize(" << FileSize << ") " << endl << flush;
    } else {
        len = pwrite64(fd,SwapData,FileSize,0L);
    }
#else
    msync(SwapData,FileSize,MS_SYNC);
#endif // WRITE_LARGE_FILE 
    munmap(SwapData,FileSize);
    close(fd);
#ifdef _REMOVE_SWAPFILE_
    if ( isDel ) unlink(SwapFileName);
#endif //_REMOVE_SWAPFILE_
    SwapData=NULL;
}


void* LoadSwapFile2D(char* FileName,
                     int x,
                     int y,
                     unsigned long DataSize,
                     int*  pf,
                     int* _fd,
                     ofstream* m_stream 
#ifdef _MPI
                     ,int rank
#endif //_MPI
                    ) {
    unsigned long FileSize=0;
    char*         data;
    void*         SwapData;
    *pf = 0;
    struct  stat  FileStat;
    int           fd;

    if ( FileName == NULL||
         x == 0          ||
         y == 0          ||
         DataSize == 0 ) {
        return(NULL);
    }

    data=new char[DataSize*x];
    memset(data,0,DataSize*x);

    SwapData=NULL;

    int yyy;
    if ( stat(FileName,&FileStat)==0 ) {
#ifdef _MPI
        if ( rank == 0 ) {
#endif //_MPI
            *m_stream << "Use old 2D-swap file \""<< FileName  << "\"..." ;
            m_stream->flush();
#ifdef _MPI
        }
#endif //_MPI
        fd = open(FileName,O_RDWR|O_SYNC,S_IRUSR|S_IWUSR);
        if ( fd == -1 ) {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Error open 2D-swap file \"" << FileName  << "\"\n" ;
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            return(NULL);
        }

        FileSize=FileStat.st_size;

        if ( FileSize!=DataSize*x*y ) {
            close(fd);
            unlink(FileName);
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "\nOld 2D-swap file has bad size.\n" ;
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            fd = open(FileName,O_CREAT|O_RDWR|O_SYNC,S_IRUSR|S_IWUSR);
            if ( fd == -1 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    *m_stream << "Error open 2D-swap file \""<< FileName  << "\"\n" ;
                    m_stream->flush();
#ifdef _MPI
                }
#endif //_MPI
                delete data;
                return(NULL);
            }

            FileSize=DataSize*x*y;

            for ( yyy=0;yyy<y;yyy++ ) {
                write(fd,data,DataSize*x);
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    *m_stream << "Create new 2D-swap file " << ((yyy+1)*x*DataSize)/1024 << " Kb ("<< 
                    ((yyy+1)*100./y) <<"%)...\r"  <<  flush;
                    m_stream->flush();
#ifdef _MPI
                }
#endif //_MPI
            }
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Create new 2D-swap file " << (y*x*DataSize)/1024 << " Kb (100%) OK.\n" ;
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
        } else
            *pf = 1;

        SwapData   = (void*)GetFileMapW(fd,FileSize
#ifdef _MPI
                                        ,rank
#endif //_MPI
                                       );
        if ( SwapData != NULL ) {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "OK\n";
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            *_fd = fd; // NEW !!!
            return(SwapData);
        } else {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Error create 2D-swap file \""<< FileName << "\"\n";
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            return(NULL);
        }
    } else {
#ifdef _MPI
        if ( rank == 0 ) {
#endif //_MPI
            *m_stream << "2D-Swap file \""<< FileName  << "\" does not exist.\nCreate new 2D-swap file\r" ;
            m_stream->flush();
#ifdef _MPI
        }
#endif //_MPI
        fd = open(FileName,O_CREAT|O_RDWR|O_SYNC,S_IRUSR|S_IWUSR);
        if ( fd == -1 ) {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Error create 2D-swap file \""<< FileName  << "\"\n" ;
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            return(NULL);
        }

        FileSize=DataSize*x*y;

        for ( yyy=0;yyy<y;yyy++ ) {
            write(fd,data,DataSize*x);
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Create new 2D-swap file " << ((yyy+1)*x*DataSize)/1024 << " Kb ("<< 
                ((yyy+1)*100./y) <<"%)...\r"  <<  flush;
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
        }
        SwapData   = (void*)GetFileMapW(fd,FileSize
#ifdef _MPI
                                        ,rank
#endif //_MPI
                                       );
        if ( SwapData != NULL ) {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Create new 2D-swap file " << (y*x*DataSize)/1024 << " Kb (100%) OK.\n" ;
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            *_fd = fd; // NEW !!!
            m_stream->flush();
            return(SwapData);
        } else {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Error create 2D-swap file \""<< FileName << "\"";
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            return(NULL);
        }
    }
}

void* LoadSwapFile1D(char* FileName,
                     int x,
                     unsigned long DataSize,
                     int*  pf,
                     int* _fd,
                     ofstream* m_stream
#ifdef _MPI
                     ,int rank
#endif //_MPI
                    ) {
    unsigned long FileSize=0;
    char*         data;
    void*         SwapData;
    *pf = 0;
    struct  stat  FileStat;
    int           fd;

    if ( FileName == NULL||
         x == 0          ||
         DataSize == 0 ) {
        return(NULL);
    }

    data=new char[DataSize];
    memset(data,0,DataSize);
    SwapData=NULL;

    int xxx;
    if ( stat(FileName,&FileStat)==0 ) {
#ifdef _MPI
        if ( rank == 0 ) {
#endif //_MPI
            *m_stream << "Use old swap file \""<< FileName  << "\"..." ;
            m_stream->flush();
#ifdef _MPI
        }
#endif //_MPI
        fd = open(FileName,O_RDWR|O_SYNC,S_IRUSR|S_IWUSR);
        if ( fd == -1 ) {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Error open swap file \"" << FileName  << "\"\n" ;
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            return(NULL);
        }

        FileSize=FileStat.st_size;

        if ( FileSize!=DataSize*x ) {
            close(fd);
            unlink(FileName);
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "\nOld swap file has bad size.\n" ;
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            fd = open(FileName,O_CREAT|O_RDWR|O_SYNC,S_IRUSR|S_IWUSR);
            if ( fd == -1 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    *m_stream << "Error open swap file \""<< FileName  << "\"\n" ;
                    m_stream->flush();
#ifdef _MPI
                }
#endif //_MPI
                delete data;
                return(NULL);
            }

            FileSize=DataSize*x;

            for ( xxx=0;xxx<x;xxx++ ) {
                write(fd,data,DataSize);
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    *m_stream << "Create new swap file " << (xxx*DataSize)/1024 << " Kb ("<< (xxx*DataSize)*100/FileSize <<"%)...\r"  ;
                    m_stream->flush();
#ifdef _MPI
                }
#endif //_MPI
            }
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Create new swap file " << (x*DataSize)/1024 << " Kb (100%) OK.\n" ;
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
        } else
            *pf = 1;

        SwapData   = (void*)GetFileMapW(fd,FileSize
#ifdef _MPI
                                        ,rank
#endif //_MPI
                                       );
        if ( SwapData != NULL ) {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "OK\n";
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            *_fd = fd; // NEW !!!
            return(SwapData);
        } else {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Error create swap file \""<< FileName << "\"\n";
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            return(NULL);
        }
    } else {
#ifdef _MPI
        if ( rank == 0 ) {
#endif //_MPI
            *m_stream << "Swap file \""<< FileName  << "\" does not exist.\nCreate new swap file\r" ;
            m_stream->flush();
#ifdef _MPI
        }
#endif //_MPI
        fd = open(FileName,O_CREAT|O_RDWR|O_SYNC,S_IRUSR|S_IWUSR);
        if ( fd == -1 ) {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Error create swap file \""<< FileName  << "\"\n" ;
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            return(NULL);
        }

        FileSize=DataSize*x;

        for ( xxx=0;xxx<x;xxx++ ) {
            write(fd,data,DataSize);
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Create new swap file " << (xxx*DataSize)/1024 << " Kb (" << (xxx*DataSize)*100/FileSize <<"%)...\r";
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
        }
        SwapData   = (void*)GetFileMapW(fd,FileSize
#ifdef _MPI
                                        ,rank
#endif //_MPI
                                       );
        if ( SwapData != NULL ) {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Create new swap file " << (x*DataSize)/1024 << " Kb (100%) OK.\n" ;
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            *_fd = fd; // NEW !!!
            m_stream->flush();
            return(SwapData);
        } else {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *m_stream << "Error create swap file \""<< FileName << "\"";
                m_stream->flush();
#ifdef _MPI
            }
#endif //_MPI
            delete data;
            return(NULL);
        }
    }
}

/* -- end file func ---*/

/* Data constructors */

Data::Data(char* name, FP val) {
    dt = DT_FLOAT;
    if ( name != 0 ) {
        Name = new char[strlen(name)+1];
        strcpy(Name,name);
    }
    fVal = val;
}

Data::Data(char* name, int val) {
    dt = DT_INT;
    if ( name != 0 ) {
        Name = new char[strlen(name)+1];
        strcpy(Name,name);
    }
    iVal = val;
}

Data::Data(char* name, char* val) {
    dt = DT_STRING;
    if ( val != NULL ) {
        StrValSize = strlen(val);
        if ( StrValSize > 0 ) {
            StrVal = new char[StrValSize+1];
            strcpy(StrVal,val);
        }
    }

    if ( name != 0 ) {
        Name = new char[strlen(name)+1];
        strcpy(Name,name);
    }

}

/* Data destructor */
Data::~Data() {
    if ( Name )
        delete [] Name;
    if ( StrValSize > 0 )
        delete [] StrVal;
    Name=NULL;
    StrVal=NULL;
}

/* Data converter */
int      Data::ConvertDataType(DataType _dt) {
    char TmpBuff[2048];
    int    i;

    if ( dt == _dt )
        return 1;
    if ( _dt == DT_STRING ) {
        if ( dt == DT_FLOAT ) {
            sprintf(TmpBuff,"%g",fVal);
            StrValSize = strlen(TmpBuff);
            StrVal = new char[StrValSize+1];
            strcpy(StrVal,TmpBuff);
            dt = DT_STRING;
        } else if ( dt == DT_INT ) {
            sprintf(TmpBuff,"%i",iVal);
            StrValSize = strlen(TmpBuff);
            StrVal = new char[StrValSize+1];
            strcpy(StrVal,TmpBuff);
            dt = DT_STRING;
        }
    } else if ( _dt == DT_FLOAT ) {
        if ( dt == DT_STRING ) {
            if ( StrVal == 0 )
                return 0;
            for ( i=0;i<StrValSize;i++ ) {
                if ( !isdigit(StrVal[i]) )
                    if ( StrVal[i]!= ' ' )
                        if ( StrVal[i]!= '-' )
                            if ( StrVal[i]!= '+' )
                                if ( StrVal[i]!= '.' )
                                    if ( StrVal[i]!= 'e' )
                                        if ( StrVal[i]!= 'E' )
                                            return 0;
            }
            fVal = atof(StrVal);
            dt = DT_FLOAT;
            StrValSize = 0;
            delete [] StrVal;
        } else if ( dt == DT_INT ) {
            fVal = (FP)(iVal);
            dt = DT_FLOAT;
        }
    } else if ( _dt == DT_INT ) {
        if ( dt == DT_STRING ) {
            if ( StrVal == 0 )
                return 0;
            for ( i=0;i<StrValSize;i++ ) {
                if ( !isdigit(StrVal[i]) )
                    if ( StrVal[i]!= ' ' )
                        if ( StrVal[i]!= '-' )
                            if ( StrVal[i]!= '+' )
                                return 0;
            }
            iVal = atoi(StrVal);
            dt = DT_INT;
            StrValSize = 0;
            delete [] StrVal;
        } else if ( dt == DT_FLOAT ) {
            iVal = (int)(fVal);
            dt = DT_INT;
        }
    }
    return 1;
}

/* Data comparition operators (< > == !=) */

int Data::operator  > (Data d) {
    if ( dt == DT_STRING )
        return strcmp(StrVal,d.StrVal);
    if ( dt == DT_FLOAT ) {
        if ( fVal > d.fVal )
            return 1;
        else
            return 0;
    } else if ( dt == DT_INT ) {
        if ( iVal > d.iVal )
            return 1;
        else
            return 0;
    }
    return 0;
}

int Data::operator  < (Data d) {
    if ( dt == DT_STRING )
        return strcmp(d.StrVal,StrVal);
    if ( dt == DT_FLOAT ) {
        if ( fVal < d.fVal )
            return 1;
        else
            return 0;
    } else if ( dt == DT_INT ) {
        if ( iVal < d.iVal )
            return 1;
        else
            return 0;
    }
    return 0;
}

int Data::operator == (Data d) {
    if ( dt == DT_STRING ) {
        if ( strcmp(StrVal,d.StrVal)==0 )
            return 1;
        else
            return 0;
    } else if ( dt == DT_FLOAT ) {
        if ( fVal == d.fVal )
            return 1;
        else
            return 0;
    } else if ( dt == DT_INT ) {
        if ( iVal == d.iVal )
            return 1;
        else
            return 0;
    }
    return 0;
}

int Data::operator != (Data d) {
    if ( dt == DT_STRING ) {
        if ( strcmp(StrVal,d.StrVal)!=0 )
            return 1;
        else
            return 0;
    } else if ( dt == DT_FLOAT ) {
        if ( fVal != d.fVal )
            return 1;
        else
            return 0;
    } else if ( dt == DT_INT ) {
        if ( iVal != d.iVal )
            return 1;
        else
            return 0;
    }
    return 0;
}


FP GetVal(XY_Table* XY, FP x ) {
    int     i, n;
    FP  y;

    //ќбщее число значений в таблице
    n = XY->n;

    //¬ таблице - единственное значение.
    if ( n == 1 )
        return( XY->y[0] );

    //јргумент меньше минимального значени€.
    if ( x <= XY->x[0] ) {
        i = 1;
        goto EndGetVal;
    }

    //јргумент больше максимального значени€.
    if ( x >= XY->x[n-1] ) {
        i = n - 1;
        goto EndGetVal;
    }

    for ( i=1; i<n; i++ ) {
        if ( (x >= XY->x[i-1]) && (x < XY->x[i]) )
            break;
    }

    EndGetVal:

    y = XY->y[i] + (XY->y[i-1] - XY->y[i])*(x - XY->x[i])/(XY->x[i-1] - XY->x[i]);

    if ( y < 0. )
        y = 0.01;

    return( y );
}

/*  InputData constructor */
InputData::InputData(char* DSource,DATA_SOURCE DS,ostream* fs,int data_size, long timeout
#ifdef _MPI
                     ,int rank
#endif //_MPI
                    ) {
    dataArray     = new UArray<Data*>();
    tableArray    = new UArray<Table*>();
    DataSource    = DSource; // File name or raw data
    DataName      = NULL;    // Data name (must store in DataSource in <start/...> directive)
    ds            = DS;
    MessageStream = fs;
    DataSize      = data_size;
    Tmp_Buff      = new char[2048];
    err_no        = 0;
    if ( ds == DS_FILE ) {
        InputDataFile = new ifstream();
        InputDataFile->open(DataSource,ios::in);
        if ( !(*InputDataFile) ) {
            err_no = -1;         
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *MessageStream << "\nFile "<< DataSource << " not found.\n";
#ifdef _MPI
            }
#endif //_MPI
        } else
            GetData(timeout
#ifdef _MPI
                    , rank
#endif //_MPI
                   );
    } else if ( ds == DS_PIPE ) {
        InputDataFile = new ifstream(0);
        GetData(timeout
#ifdef _MPI
                , rank
#endif //_MPI
               );
    } else if ( ds == DS_MEM )
        GetData(timeout
#ifdef _MPI
                , rank
#endif //_MPI
               );
}

/*  GetData */
int InputData::GetData(long timeout
#ifdef _MPI
                       ,int rank
#endif //_MPI
                      ) {
    int ret = 0;

    if ( ds == DS_FILE )
        ret = GetDataFromFile(timeout
#ifdef _MPI
                              , rank
#endif //_MPI
                             );
    else if ( ds == DS_PIPE )
        ret = GetDataFromFile(timeout
#ifdef _MPI
                              , rank
#endif //_MPI
                             );
    else if ( ds == DS_MEM )
        ret = GetDataFromMemory(
#ifdef _MPI
                               rank
#endif //_MPI
                               );
    err_no = ret;
    return err_no;
}


int  InputData::GetDataFromMemory(
#ifdef _MPI
                                 int rank
#endif //_MPI
                                 ) {
    char*      TmpBuff=NULL;
    char*      EndBuff=NULL;
    char*      TmpPtr =NULL;
    char*      TmpPtr1=NULL;
    char       LocalName[2048];  
    unsigned int buff_size,i,s=0;
    char*      GlobalPtr=DataSource;               
    char*      PtrFlag=DataSource;

    Data*      D;
    Table*     T;

    if ( DataSize != 0 )
        buff_size = DataSize;
    else
        buff_size = strlen(DataSource);

    TmpBuff   = new char[2048];
#ifdef _MPI
    if ( rank == 0 ) {
#endif //_MPI
        *MessageStream << " \n" ;
        MessageStream->flush();
#ifdef _MPI
    }
#endif //_MPI

    do {
        memset(TmpBuff,0,1024);

        PtrFlag=(char*)memccpy(TmpBuff,GlobalPtr,'\n',min(1024,buff_size));
        buff_size -= strlen(TmpBuff);
        GlobalPtr += strlen(TmpBuff);
#ifdef _DEBUG_100
#ifdef _MPI
        if ( rank == 0 ) {
#endif //_MPI
            *MessageStream << TmpBuff ;
            MessageStream->flush();
#ifdef _MPI
        }
#endif //_MPI
#endif
        strtok(TmpBuff,"#;");
        TmpPtr  = strstr(TmpBuff,"<start/");
        if ( TmpPtr != 0 )
            TmpPtr1 = strtok(TmpPtr,">");
        if ( TmpPtr != 0 ) {
            if ( s==1 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\n<start/...> directive define twice:\n%s\n",TmpBuff);
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                if ( TmpBuff !=0 )
                    delete  TmpBuff;
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 2;
            }
            DataName = new char[strlen(TmpPtr+7)+1];
            strcpy(DataName,TmpPtr+7);
            EndBuff = new char[strlen(DataName)+20];
            sprintf(EndBuff,"<end/%s>",DataName);
            sprintf(Message,"Load \"%s\" data...",DataName);
            *MessageStream << Message ;
            MessageStream->flush();
            s++;
        }
        TmpPtr = strstr(TmpBuff,"<data/");

        if ( TmpPtr != 0 ) {
            if ( s==0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\n<start/...> directive not found.\n");
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                if ( TmpBuff !=0 )
                    delete  TmpBuff;
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 2;
            }
            TmpPtr1 = strstr(TmpPtr+6,"=");
            if ( TmpPtr1 == 0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\nError <data/...> directive:\n %s\n",TmpBuff);
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                if ( TmpBuff !=0 )
                    delete  TmpBuff;
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 3;
            }
            TmpPtr1 = strtok(TmpPtr+6,"=");
            strcpy(LocalName,TmpPtr+6);
            TmpPtr1 = strtok(TmpPtr+7+strlen(LocalName),">");
            if ( TmpPtr1 == 0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\nError <data/...> directive:\n %s\n",TmpBuff);
#ifdef _MPI
                }
#endif //_MPI
                *MessageStream << Message ;
                MessageStream->flush();
                if ( TmpBuff !=0 )
                    delete  TmpBuff;
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 3;
            }
            D = new Data(LocalName,TmpPtr1);
            dataArray->AddElement(&D);
        }
        TmpPtr = strstr(TmpBuff,"<table=");
        if ( TmpPtr != 0 ) {
            if ( s==0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\n<start/...> directive not found.\n");
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                if ( TmpBuff !=0 )
                    delete  TmpBuff;
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 2;
            }
            TmpPtr1 = strstr(TmpPtr+7,"/");
            if ( TmpPtr1 == 0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\nError <table=.../...> directive:\n %s\n",TmpBuff);
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                if ( TmpBuff !=0 )
                    delete  TmpBuff;
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 3;
            }
            TmpPtr1 = strstr(TmpPtr,"=");
            if ( TmpPtr != 0 ) {
                strtok(TmpPtr1,"/");
                strcpy(LocalName,TmpPtr1+1);
                TmpPtr1=strtok(0,">");
            } else if ( TmpPtr1 == 0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\nError <table=.../...> directive:\n %s\n",TmpBuff);
                    *MessageStream << Message ;
                    MessageStream->flush(); 
#ifdef _MPI
                }
#endif //_MPI

                if ( TmpBuff !=0 )
                    delete  TmpBuff;
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 3;
            }

            T = new Table(LocalName,atoi(TmpPtr1));

            for ( i=0; i<T->n; i++ ) {
                memset(TmpBuff,0,1024);

                PtrFlag=(char*)memccpy(TmpBuff,GlobalPtr,'\n',min(1024,buff_size));

#ifdef _DEBUG_100
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    *MessageStream  << TmpBuff << "\n" ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
#endif //_DEBUG_100 

                buff_size -= strlen(TmpBuff);
                GlobalPtr += strlen(TmpBuff);

                if ( strstr(TmpBuff,"<endtable>")!=0 )
                    break;
                T->x[i] = atof(TmpBuff);
                TmpPtr1 = strstr(TmpBuff,"\t");
                if ( TmpPtr1 == 0 )
                    TmpPtr1 = strstr(TmpBuff," ");
                if ( TmpPtr1 == 0 ) {
#ifdef _MPI
                    if ( rank == 0 ) {
#endif //_MPI
                        sprintf(Message,"\nError <table=.../...> directive:\n %s\n",TmpBuff);
                        *MessageStream << Message ;
                        MessageStream->flush();
#ifdef _MPI
                    }
#endif //_MPI
                    if ( TmpBuff !=0 )
                        delete  TmpBuff;
                    if ( EndBuff !=0 )
                        delete  EndBuff;
                    return 3;
                }
                T->y[i] = atof(TmpPtr1);
            }
            tableArray->AddElement(&T);
        }

        if ( EndBuff !=0 )
            if ( TmpBuff !=0 )
                if ( strstr(TmpBuff,EndBuff)!=0 ) {
                    if ( TmpBuff !=0 )
                        delete  TmpBuff;
                    if ( EndBuff !=0 )
                        delete  EndBuff;
#ifdef _MPI
                    if ( rank == 0 ) {
#endif //_MPI
                        sprintf(Message,"OK\n");
                        *MessageStream << Message ;
                        MessageStream->flush(); 
#ifdef _MPI
                    }
#endif //_MPI
                    return 0;
                }
    }while ( PtrFlag );

    if ( TmpBuff !=0 )
        delete  TmpBuff;
    if ( EndBuff !=0 )
        delete  EndBuff;
    if ( s==0 ) {
#ifdef _MPI
        if ( rank == 0 ) {
#endif //_MPI
            sprintf(Message,"\n<start/...> directive not found.\n");
            *MessageStream << Message ;
            MessageStream->flush();
#ifdef _MPI
        }
#endif //_MPI
        return 2;
    }
#ifdef _MPI
    if ( rank == 0 ) {
#endif //_MPI
        sprintf(Message,"\n<end/...> directive not found.\n");
        *MessageStream << Message ;
        MessageStream->flush();
#ifdef _MPI
    }
#endif //_MPI
    return 4;
}

int  InputData::SaveAllDataAsText(char* OutBuff) {
    return strlen(OutBuff);
}


int  InputData::GetDataFromFile(long timeout
#ifdef _MPI
                                ,int rank
#endif //_MPI
                               ) {
    char*      EndBuff=NULL;
    char*      TmpPtr =NULL;
    char*      TmpPtr1=NULL;
    char       LocalName[2048];
    unsigned int  i,s=0;
    Data*      D;
    Table*     T;
    fd_set     rfds;
    struct     timeval tv;
    int        retval;
    if ( !InputDataFile->is_open() ) {
#ifdef _MPI
        if ( rank == 0 ) {
#endif //_MPI
            *MessageStream << "Error open data file \"" << DataSource << "\".\n" ;
            MessageStream->flush();
#ifdef _MPI
        }
#endif //_MPI
        return 1;
    }
    do {
        memset(Tmp_Buff,0,1024);
        if ( ds==DS_PIPE ) {
/* Watch InputDataFile->filedesc() to see when it has input. */
            tv.tv_sec = timeout;
            tv.tv_usec = 0;
            FD_ZERO(&rfds);
            FD_SET(fd, &rfds); 
/* Wait up to *timeout* seconds (default timeout=1 sec). */
            retval = select(fd+1, &rfds, NULL, NULL, &tv);
/* Don't rely on the value of tv now! */
        } else {
            retval = 1;
        }

        if ( retval )
            InputDataFile->getline(Tmp_Buff,1000,'\n');
        else {
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                *MessageStream << "Error read data steram, timeout expired (" << timeout << " sec).\n";
                MessageStream->flush();
#ifdef _MPI
            }
#endif //_MPI
            return 4;
        }
#ifdef _DEBUG_100
#ifdef _MPI
        if ( rank == 0 ) {
#endif //_MPI
            sprintf(Message,"%s",Tmp_Buff);
            *MessageStream << Message ;
            MessageStream->flush();
#ifdef _MPI
        }
#endif //_MPI
#endif //_DEBUG_100 

        strtok(Tmp_Buff,"#;");
        TmpPtr  = strstr(Tmp_Buff,"<start/");
        if ( TmpPtr != 0 )
            TmpPtr1 = strtok(TmpPtr,">");
        if ( TmpPtr != 0 ) {
            if ( s==1 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\n<start/...> directive define twice:\n%s\n",Tmp_Buff);
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                InputDataFile->close();
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 2;
            }
            DataName = new char[strlen(TmpPtr+7)+1];
            strcpy(DataName,TmpPtr+7);
            EndBuff = new char[strlen(DataName)+20];
            sprintf(EndBuff,"<end/%s>",DataName);
#ifdef _MPI
            if ( rank == 0 ) {
#endif //_MPI
                sprintf(Message,"Load \"%s\" data...",DataName);
                *MessageStream << Message ;
                MessageStream->flush();
#ifdef _MPI
            }
#endif //_MPI
            s++;
        }
        TmpPtr = strstr(Tmp_Buff,"<data/");

        if ( TmpPtr != 0 ) {
            if ( s==0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\n<start/...> directive not found.\n");
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                InputDataFile->close();
                //if(Tmp_Buff !=0 )
                //   delete  Tmp_Buff;
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 2;
            }
            TmpPtr1 = strstr(TmpPtr+6,"=");
            if ( TmpPtr1 == 0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\nError <data/...> directive:\n %s\n",Tmp_Buff);
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                InputDataFile->close();
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 3;
            }
            TmpPtr1 = strtok(TmpPtr+6,"=");
            strcpy(LocalName,TmpPtr+6);
            TmpPtr1 = strtok(TmpPtr+7+strlen(LocalName),">");
            if ( TmpPtr1 == 0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\nError <data/...> directive:\n %s\n",Tmp_Buff);
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                InputDataFile->close();
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 3;
            }
            D = new Data(LocalName,TmpPtr1);
            dataArray->AddElement(&D);
        }
        TmpPtr = strstr(Tmp_Buff,"<table=");
        if ( TmpPtr != 0 ) {
            if ( s==0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\n<start/...> directive not found.\n");
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                InputDataFile->close();
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 2;
            }
            TmpPtr1 = strstr(TmpPtr+7,"/");
            if ( TmpPtr1 == 0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\nError <table=.../...> directive:\n %s\n",Tmp_Buff);
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                InputDataFile->close();
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 3;
            }
            TmpPtr1 = strstr(TmpPtr,"=");
            if ( TmpPtr != 0 ) {
                strtok(TmpPtr1,"/");
                strcpy(LocalName,TmpPtr1+1);
                TmpPtr1=strtok(0,">");
            } else if ( TmpPtr1 == 0 ) {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    sprintf(Message,"\nError <table=.../...> directive:\n %s\n",Tmp_Buff);
                    *MessageStream << Message ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                InputDataFile->close();
                if ( EndBuff !=0 )
                    delete  EndBuff;
                return 3;
            }
            T = new Table(LocalName,atoi(TmpPtr1));

            for ( i=0; i<T->n; i++ ) {
                memset(Tmp_Buff,0,1024);
                if ( ds==DS_PIPE ) {
/* Watch InputDataFile->filedesc() to see when it has input. */
                    FD_ZERO(&rfds);
/* Wait up to *timeout* seconds (default timeout=1 sec). */
                    tv.tv_sec = timeout;
                    tv.tv_usec = 0;
                    FD_SET(fd, &rfds); // For icc only !!!
                    retval = select(fd+1, &rfds, NULL, NULL, &tv);
/* Don't rely on the value of tv now! */
                } else {
                    retval = 1;
                }

                if ( retval )
                    InputDataFile->getline(Tmp_Buff,1024);
                else
                    return 4;

                if ( strstr(Tmp_Buff,"<endtable>")!=0 )
                    break;
                T->x[i] = atof(Tmp_Buff);
                TmpPtr1=strstr(Tmp_Buff," ");
                if ( TmpPtr1 == 0 ) {
#ifdef _MPI
                    if ( rank == 0 ) {
#endif //_MPI
                        sprintf(Message,"\nError <table=.../...> directive:\n %s\n",Tmp_Buff);
                        *MessageStream << Message ;
                        MessageStream->flush();
#ifdef _MPI
                    }
#endif //_MPI
                    InputDataFile->close();
                    //if(Tmp_Buff !=0 )
                    //   delete  Tmp_Buff;
                    if ( EndBuff !=0 )
                        delete  EndBuff;
                    return 3;
                }
                T->y[i] = atof(TmpPtr1);
            }
            tableArray->AddElement(&T);
        }

        if ( EndBuff !=0 )
            if ( Tmp_Buff !=0 )
                if ( strstr(Tmp_Buff,EndBuff)!=0 ) {
                    //if(Tmp_Buff !=0 )
                    //   delete  Tmp_Buff;
                    if ( EndBuff !=0 )
                        delete  EndBuff;
                    InputDataFile->close();
#ifdef _MPI
                    if ( rank == 0 ) {
#endif //_MPI

                        sprintf(Message,"OK\n");
                        *MessageStream << Message ;
                        MessageStream->flush(); 
#ifdef _MPI
                    }
#endif //_MPI

                    return 0;
                }
    }while ( InputDataFile );

    //if(Tmp_Buff !=0 )
    //   delete  Tmp_Buff;
    if ( EndBuff !=0 )
        delete  EndBuff;
    if ( s==0 ) {
#ifdef _MPI
        if ( rank == 0 ) {
#endif //_MPI
            sprintf(Message,"<start/...> directive not found.\n");
            *MessageStream << Message ;
            MessageStream->flush();
#ifdef _MPI
        }
#endif //_MPI
        InputDataFile->close();
        return 2;
    }
#ifdef _MPI
    if ( rank == 0 ) {
#endif //_MPI
        sprintf(Message,"\n<end/...> directive not found.\n");
        *MessageStream << Message ;
        MessageStream->flush();
#ifdef _MPI
    }
#endif //_MPI
    //InputDataFile->close();
    //delete InputDataFile;
    return 4;
}

void InputData::ClearInputData() {
    unsigned int i;

    if ( DataName )
        delete DataName;
    else
        return;

    unsigned int id = dataArray->GetNumElements();
    unsigned int it = tableArray->GetNumElements();

    for ( i=0;i<id;i++ )
        delete dataArray->GetElement(i);

    dataArray->CleanArray();

    for ( i=0;i<it;i++ )
        delete  tableArray->GetElement(i);

    tableArray->CleanArray();

    DataName=NULL;
}

InputData::~InputData() {
    ClearInputData();
    delete dataArray;
    delete tableArray;
    delete Tmp_Buff;
    if ( ds == DS_FILE )
        InputDataFile->close();
    delete InputDataFile;
}

ostream* InputData::GetMessageStream() {
    return MessageStream;
}
void      InputData::SetMessageStream(ostream* fs) {
    MessageStream=fs;
}

DATA_SOURCE InputData::GetDataSource() {
    return ds;
}
int         InputData::GetDataError() {
    return err_no;
}
char*       InputData::GetDataName() {
    return DataName;
}

int   InputData::GetIntVal(char* Name
#ifdef _MPI
                           ,int rank
#endif //_MPI
                          ) {
    int      i,retVal;
    Data*    D;
    DataType SaveDT;

    for ( i=0;i<(int)dataArray->GetNumElements();i++ ) {
        D = *(dataArray->GetElementPtr(i));
        if ( strcmp(D->GetName(),Name)==0 ) {
            SaveDT = D->GetDataType();
            D->ConvertDataType(DT_INT);

            if ( D->GetDataType()==DT_INT ) {
                retVal =  D->iVal;
                D->ConvertDataType(SaveDT);
                return retVal;
            } else {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    *MessageStream << "Data object \"" << Name << "\" have not "<< _INT_TYPE<<" type.\n" ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                err_no = -1;
                return 0;
            }
        }
    }
#ifdef _MPI
    if ( rank == 0 ) {
#endif //_MPI
        *MessageStream << "Data object \"" << Name << "\" not found in \"" << InputData::DataName << "\".\n";
        MessageStream->flush(); 
        *MessageStream << "Please, check input data file.\n" ;
        MessageStream->flush();
#ifdef _MPI
    }
#endif //_MPI
    err_no = -1;
    return 0;
}

FP    InputData::GetFloatVal(char* Name
#ifdef _MPI
                                 ,int rank
#endif //_MPI
                                ) {
    int    i;
    Data*  D;
    DataType SaveDT;
    FP retVal;

    for ( i=0;i<(int)dataArray->GetNumElements();i++ ) {
        D = *(dataArray->GetElementPtr(i));
        if ( strcmp(D->GetName(),Name)==0 ) {
            SaveDT = D->GetDataType();
            D->ConvertDataType(DT_FLOAT);

            if ( D->GetDataType()==DT_FLOAT ) {
                retVal =  D->fVal;
                D->ConvertDataType(SaveDT);
                return retVal;
            } else {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    *MessageStream << "Data object \"" << Name << "\" have not "<< _FLOAT_TYPE<<" type.\n" ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                err_no = -1;
                return  0.;
            }
        }
    }
#ifdef _MPI
    if ( rank == 0 ) {
#endif //_MPI
        *MessageStream << "Data object \"" << Name << "\" not found in \"" << InputData::DataName << "\".\n";
        *MessageStream << "Please, check input data file.\n" ;
        MessageStream->flush();
#ifdef _MPI
    }
#endif //_MPI
    err_no = -1;
    return 0.;
}

char*     InputData::GetStringVal(char* Name
#ifdef _MPI
                                  ,int rank
#endif //_MPI
                                 ) {
    int    i;
    Data*  D;
    DataType SaveDT;

    for ( i=0;i<(int)dataArray->GetNumElements();i++ ) {
        D = *(dataArray->GetElementPtr(i));
        if ( strcmp(D->GetName(),Name)==0 ) {
            SaveDT = D->GetDataType();
            D->ConvertDataType(DT_STRING);

            if ( D->GetDataType()==DT_STRING ) {
                return D->StrVal;
            } else {
#ifdef _MPI
                if ( rank == 0 ) {
#endif //_MPI
                    *MessageStream << "Data object \"" << Name << "\" have not "<< _STRING_TYPE<<" type.\n" ;
                    MessageStream->flush();
#ifdef _MPI
                }
#endif //_MPI
                err_no = -1;
                return 0;
            }
        }
    }
#ifdef _MPI
    if ( rank == 0 ) {
#endif //_MPI
        *MessageStream << "Data object \"" << Name << "\" not found in \"" << InputData::DataName << "\".\n";
        *MessageStream << "Please, check input data file.\n" ;
        MessageStream->flush();
#ifdef _MPI
    }
#endif //_MPI
    err_no = -1;
    return 0;
}

DataType  InputData::GetDataType(char* Name
#ifdef _MPI
                                 ,int rank
#endif //_MPI
                                ) {
    int    i;
    Data*  D;

    for ( i=0;i<(int)dataArray->GetNumElements();i++ ) {
        D = *(dataArray->GetElementPtr(i));
        if ( strcmp(D->GetName(),Name)==0 )
            return D->GetDataType();
    }
#ifdef _MPI
    if ( rank == 0 ) {
#endif //_MPI
        *MessageStream << "Data object \"" << Name << "\" not found in \"" << InputData::DataName << "\".\n";
        *MessageStream << "Please, check input data file.\n" ;
        MessageStream->flush();
#ifdef _MPI
    }
#endif //_MPI
    err_no = -1;
    return DT_UNKNOWN;
}


FP  InputData::GetVal(char* Name, FP par
#ifdef _MPI
                          ,int rank
#endif //_MPI
                         ) {
    int    i;
    Table* T;

    for ( i=0;i<(int)tableArray->GetNumElements();i++ ) {
        T = *(tableArray->GetElementPtr(i));
        if ( strcmp(T->GetName(),Name)==0 )
            return T->GetVal((FP)par);
    }
#ifdef _MPI
    if ( rank == 0 ) {
#endif //_MPI
        *MessageStream << "Data object \"" << Name << "\" not found in \"" << InputData::DataName << "\".\n";
        *MessageStream << "Please, check input data file.\n" ;
        MessageStream->flush();
#ifdef _MPI
    }
#endif //_MPI
    err_no = -1;
    return 0.;
}

Table* InputData::zeroTable = new Table((char*)"ZeroTable",0);

Table*  InputData::GetZeroTable() {
    return zeroTable;
}

Table*  InputData::GetTable(char* Name
#ifdef _MPI
                            ,int rank
#endif //_MPI
                           ) {
    int    i;
    Table* T;

    for ( i=0;i<(int)tableArray->GetNumElements();i++ ) {
        T = *(tableArray->GetElementPtr(i));
        if ( strcmp(T->GetName(),Name)==0 )
            return T;
    }
#ifdef _MPI
    if ( rank == 0 ) {
#endif //_MPI
        *MessageStream << "Table object \"" << Name << "\" not found in \"" << InputData::DataName << "\".\n";
        *MessageStream << "Please, check input data file.\n" ;
        MessageStream->flush();
#ifdef _MPI
    }
#endif //_MPI
    err_no = -1;
    return zeroTable;
}

int     InputData::GetTableSize(char* Name
#ifdef _MPI
                                ,int rank
#endif //_MPI
                               ) {
    int    i;
    Table* T;

    for ( i=0;i<(int)tableArray->GetNumElements();i++ ) {
        T = *(tableArray->GetElementPtr(i));
        if ( strcmp(T->GetName(),Name)==0 )
            return T->n;
    }
#ifdef _MPI
    if ( rank == 0 ) {
#endif //_MPI
        *MessageStream << "Table object \"" << Name << "\" not found in \"" << InputData::DataName << "\".\n";
        *MessageStream << "Please, check input data file.\n" ;
        MessageStream->flush();
#ifdef _MPI
    }
#endif //_MPI
    err_no = -1;
    return -1;
}

int     InputData::EnumData(UArray<char*>* tD) {
    int    i;
    Data*  D;
    char*  td;

    for ( i=0;i<(int)dataArray->GetNumElements();i++ ) {
        D = *(dataArray->GetElementPtr(i));
        td=D->GetName();
        tD->AddElement(&td);
    }
    return tD->GetNumElements();
}

int     InputData::EnumTable(UArray<char*>* tT) {
    int    i;
    Table* T;
    char*  tt;

    for ( i=0;i<(int)tableArray->GetNumElements();i++ ) {
        T = *(tableArray->GetElementPtr(i));
        tt = T->GetName();
        tT->AddElement(&tt);
    }
    return tT->GetNumElements();
}

XY_Table::XY_Table(unsigned int N) {
    n = N;
    x = new FP[n];
    y = new FP[n];
}

XY_Table::~XY_Table() {
    if ( x != NULL )
        delete [] x;
    if ( y != NULL )
        delete [] y;
}


FP XY_Table::GetX(unsigned  int i) {
    if ( i >=n ) return 0;
    return x[i];
}

FP XY_Table::GetY(unsigned int i) {
    if ( i >=n ) return 0;
    return y[i];
}

unsigned int XY_Table::GetNumNodes() {
    return n;
}

Table::Table(char* name ,int N):XY_Table(N) {
    n = N;
    Name=NULL;
    if ( name != 0 && strlen(name)>0 ) {
        Name = new char[strlen(name)+1];
        strcpy(Name,name);
    }
}

Table::~Table() {
    if ( Name )
        delete[] Name;
    Name=NULL;
}

char*  Table::GetName() {
    return Name;
}

void   Table::SetName(char* name) {
    if ( Name ) delete[] Name;
    Name = new char[strlen(name)+1];
    if ( name ) strcpy(Name,name);
}

int    Table::operator  > (Table) {
    return 0;
}
int    Table::operator  < (Table) {
    return 0;
}

FP Table::GetVal(FP _x ) {
    if ( this == InputData::GetZeroTable() )
        return 0.;
    //ќбщее число значений в таблице
    register int    i, _n = n;
    register FP _y;

    //¬ таблице - единственное значение.
    if ( _n == 1 )
        return( y[0] );

    //јргумент меньше минимального значени€.
    if ( _x <= x[0] ) {
        i = 1;
        goto EndGetVal;

    }

    //јргумент больше максимального значени€.
    if ( _x >= x[n-1] ) {
        i = _n - 1;
        goto EndGetVal;
    }

    for ( i=1; i<_n; i++ ) {
        if ( (_x >= x[i-1]) && (_x < x[i]) )
            break;
    }

    EndGetVal:

    _y = y[i] + (y[i-1] - y[i])*(_x - x[i])/(x[i-1] - x[i]);

    /*
    if ( _y < 0. ) // ??
        _y = 0.01;
    */
    return( _y );
}


