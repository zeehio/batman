// written by Dr. Jie Hao, Dr William Astle
#include "myheader.h"


void read_txtf(matrix *data, char filename[])
{
	int BUFF_SIZE=40000;
	char buffer[BUFF_SIZE];
    char *p; // p is the pointer point to the first character of buffer

	int j=0;// i and j are row and column indeces of c, which are start from 0 to 2503 and 99, respectively
	int count=0; // count for the number of ' '
    int col = 0;
    FILE *fp=fopen(filename,"r");

    if( !fp)
    {
    printf("Can't open file %s, exiting ... \n", filename);
    system("PAUSE");
    exit(1);
    }


    char *tp0 = fgets(buffer, BUFF_SIZE,fp);
    p = buffer;
    while (*p!='\n')
    {
     p++;
     if (*p == '\t')
     col++;
    }
    

    (*data).resize(col+1);
    //resize the buffer size;
    BUFF_SIZE = col*20;
    printf("Buffer size is adjusted to %i \n", BUFF_SIZE);
    
    
    while( 1 )
    {
    char buffer[BUFF_SIZE];
      //printf("buffer created successfully. \n");
         char buffer1[BUFF_SIZE];
         //printf("buffer1 created successfully. \n");
         char * tp0 = fgets(buffer, BUFF_SIZE, fp);
      p = buffer; // the pointer 'p' point to the address of the first character of buffer
      if(feof(fp))
      break;
      while (*p != '\n')
      {
            if(*p == '\t')// space or not?
    		{
            buffer1[j]='\0';
            (*data)[count].push_back(atof(buffer1)); // convert string to float
   			count++;
   			j = 0;
   			p++;
    		}
            else
            {
            buffer1[j]= *p; // to be stored in column 1
	        j++;
      	    p++;
            }
    }
    if(*p == '\n')
    	    {
            buffer1[j] = '\0';
            (*data)[count].push_back(atof(buffer1));
            j = 0;
    	    }
    count = 0;
    j=0;
    }
    fclose(fp);

}

