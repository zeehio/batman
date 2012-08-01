#include "myheader.h"

int read_datf(vector<string> * name, matrix *c1, matrix *c2, matrix *c3, matrix *c4, 
matrix *c5, matrix *c6, matrix *c7, matrixI *c8, char filename[])
{
	char *p; // p is the pointer point to the first character of buffer 

	int i=0, j=0; // i and j are row and column indeces of c, which are start from 0 to 2503 and 99, respectively
	int count=0; // count for the number of '\t'

    FILE *fp=fopen(filename,"r");
    
    if( !fp)
    {
        printf("Can't open file %s, exiting ...\n", filename);
        system("PAUSE");
        exit(1);
    }
    //fgets(buffer, 3000, fp);
    while( 1 )
    {
    char buffer[3000] = {'\0'};
    char buffer1[100]={'\0'}; // 
    char buffer2[100]={'\0'}; // 
    char buffer3[100]={'\0'}; // 
    char buffer4[100]={'\0'}; // 
    char buffer5[100]={'\0'}; // 
    char buffer6[100]={'\0'}; // 
    char buffer7[100]={'\0'}; // 
    char buffer8[30]={'\0'}; //
    //vector<int> m(5,0);
      fgets(buffer, 3000, fp);
      p = buffer; // the pointer 'p' point to the address of the first character of buffer
      if(feof(fp))
      break;
      while( *p != '\n')  
      {
	  if (*p=='"')
	  {p++;	
	  }
      	  else if(*p == '\t')// tab or not?
    		{
    			switch (count)
    			{	
    			case 0:
    				(*name)[i]+='\0'; // end of the string
    				break;
    			case 1:
    				buffer1[j]='\0';
    				(*c1)[i].push_back(atof(buffer1)); // convert string to float
    				break;
    			case 2:
    				buffer2[j]='\0';
    				(*c2)[i].push_back(atof(buffer2));
    				break;
    			case 3:
    				buffer3[j]='\0';
    				(*c3)[i].push_back(atof(buffer3));
    				break;
    			case 4:
    				buffer4[j]='\0';
    				(*c4)[i].push_back(atof(buffer4));
    				break;
    		    case 5:
    				buffer5[j]='\0';
    				(*c5)[i].push_back(atof(buffer5));
    				break;
    		    case 6:
    				buffer6[j]='\0';
    				(*c6)[i].push_back(atof(buffer6));
    				break;
    			case 7:
    				buffer7[j]='\0';
    				(*c7)[i].push_back(atof(buffer7));
    				break;
    			}
    			j=0;
    			count++;
    			p++;
    		}
    	   else if(*p == ',' && count != 0)// , or not?
    		{
    			switch (count)
        	    {
        		case 1:
        			buffer1[j]='\0';
                    (*c1)[i].push_back(atof(buffer1));
                    break;
        		case 2:
                    buffer2[j]='\0';
        			(*c2)[i].push_back(atof(buffer2)); // column 2
        			break;
        		case 3:
                    buffer3[j]='\0';
        			(*c3)[i].push_back(atof(buffer3)); // column 3
        			break;
        		case 4:
                    buffer4[j]='\0';
        			(*c4)[i].push_back(atof(buffer4)); // to be stored in column 4
                	break;
        		case 5:
                    buffer5[j]='\0';
        			(*c5)[i].push_back(atof(buffer5)); 
        			break;
        		case 6:
                    buffer6[j]='\0';
        			(*c6)[i].push_back(atof(buffer6)); 
        			break;
        		case 7:
                    buffer7[j]='\0';
        			(*c7)[i].push_back(atof(buffer7)); 
                	break;
                case 8:
                    buffer8[j]='\0';
        			(*c8)[i].push_back(atoi(buffer8)); 
                	break;
        	    }
        	  j = 0;
    		  p++;
    		}
    	  else
            {
        	  switch (count)
        	  {
        		case 0:
        			(*name)[i]+= *p; // name 
        			break;
        		case 1:
        			buffer1[j]= *p; // to be stored in column 1
        			break;
        		case 2:
        			buffer2[j] = *p; // column 2
        			break;
        		case 3:
        			buffer3[j] = *p; // column 3
        			break;
        		case 4:
        			buffer4[j] = *p;
        			break;
        		case 5:
        			buffer5[j] = *p; 
        			break;
        		case 6:
        			buffer6[j] = *p; 
        			break;
        		case 7:
        			buffer7[j] = *p;
        			break;
        		case 8:
        			buffer8[j] = *p;
        			break;
        	  }
        	  j++;
        	  p++;
           }
        }
        if(*p == '\n')
    	 {		
            buffer8[j] = '\0';
            (*c8)[i].push_back(atoi(buffer8)); // convert strinig to int
    	 } 
      i++;
      count = 0;
      j=0;
    }
    fclose(fp);
    return i;
}
