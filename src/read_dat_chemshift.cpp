#include "myheader.h"

int read_dat_chemshift(vector<string> * name, matrix *c1,matrix *c2, char filename[], int s)
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
		char buffer[10000] = {'\0'};
		char buffer1[100]={'\0'}; // 
		char buffer2[100]={'\0'}; // 
		//char buffer3[100]={'\0'}; // 
		//char buffer4[100]={'\0'}; // 
		//char buffer5[100]={'\0'}; // 
		//char buffer6[100]={'\0'}; // 
		//char buffer7[100]={'\0'}; // 
		//char buffer8[30]={'\0'}; //
		//vector<int> m(5,0);
		fgets(buffer, 10000, fp);
		p = buffer; // the pointer 'p' point to the address of the first character of buffer
		if(feof(fp))
			break;
		while( *p != '\n')  
		{
			if (*p=='"')
			{
				p++;	
			}
			else if(*p == '\t')// tab or not?
			{
				if (count == 0)
				{	
					(*name)[i]+='\0'; // end of the string
				}	else if (count == 1) {
					buffer1[j]='\0';
					(*c1)[i].push_back(atof(buffer1));
				} else if (count == (s+2)) {
					buffer2[j]='\0';
					(*c2)[i].push_back(atof(buffer2));
				}
				j=0;
				count++;
				p++;
			}
			else
			{
				if (count == 0)
				{
					(*name)[i]+= *p; // name 
				} else if (count == 1) {
					buffer1[j] = *p;
				} else if (count == (s+2)) {
					buffer2[j] = *p;
				}
				j++;
				p++;
			}
		}
		if(*p == '\n' && count == (s+2))
		{		
			buffer2[j] = '\0';
			(*c2)[i].push_back(atof(buffer2)); // convert strinig to double at end of line
		} 
		i++;
		count = 0;
		j=0;
	}
	fclose(fp);
	return i;
}
