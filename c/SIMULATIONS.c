#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

void configparser(char configfile[20]);

int main(void){
    configparser("config.ini");
    return 0;
}
    
void configparser(char configfile[20])
{
    FILE* fp = fopen(configfile, "r");
    char buf[150], tmpbuf[150], filename[10], lightcurvetype[10];
    float nsources = 0, flmin = 0, flmax = 0, flerr = 0, dmin = 0, dmax = 0, detthresh = 0, extrathresh = 0, obssens = 0, obssig = 0, obsinterval = 0, obsduration = 0;
    int nobs = 0, countF = 0, countS = 0, countI = 0;
    
    while (fgets(buf, sizeof buf, fp) != NULL){
        if (buf[0] == 'F'){
            for(int i=1; i<150; i++){
                if((buf[i] =='\n') || (buf[i] == ' ')) break;
                tmpbuf[i-1] = buf[i];
            }
            switch(countF){
                case 0 : nsources = atof(tmpbuf);
                memset(tmpbuf,0,sizeof(tmpbuf));
                countF +=1;
                break;
                case 1 : flmin = atof(tmpbuf);
                memset(tmpbuf,0,sizeof(tmpbuf));
                countF +=1;
                break;
                case 2 : flmax = atof(tmpbuf);
                memset(tmpbuf,0,sizeof(tmpbuf));
                countF +=1;
                break;
                case 3 : flerr = atof(tmpbuf);
                memset(tmpbuf,0,sizeof(tmpbuf));
                countF +=1;
                break;
                case 4 : dmin = atof(tmpbuf);
                memset(tmpbuf,0,sizeof(tmpbuf));
                countF +=1;
                break;
                case 5 : dmax = atof(tmpbuf);
                memset(tmpbuf,0,sizeof(tmpbuf));
                countF +=1;
                break;
                case 6 : detthresh = atof(tmpbuf);
                memset(tmpbuf,0,sizeof(tmpbuf));
                countF +=1;
                break;
                case 8 : obssens = atof(tmpbuf);
                memset(tmpbuf,0,sizeof(tmpbuf));
                countF +=1;
                break;
                case 9 : obssig = atof(tmpbuf);
                memset(tmpbuf,0,sizeof(tmpbuf));
                countF +=1;
                break;
                case 10 : obsinterval = atof(tmpbuf);
                memset(tmpbuf,0,sizeof(tmpbuf));
                countF +=1;
                break;
                case 11 : obsduration = atof(tmpbuf);
                memset(tmpbuf,0,sizeof(tmpbuf));
                countF +=1;
                break;
                }
        }
        if (buf[0] == 'S'){
            for(int i=1; i<150; i++){
                if((buf[i] =='\n') || (buf[i] == ' ')) break;
                tmpbuf[i-1] = buf[i];
            }
            switch(countS){
                case 0 : for(int i=0; i<10; i++) filename[i] = tmpbuf[i];
                memset(tmpbuf,0,sizeof(tmpbuf));
                countS +=1;
                break;
                case 1 : for(int i=0; i<10; i++) lightcurvetype[i] = tmpbuf[i];
                memset(tmpbuf,0,sizeof(tmpbuf));
                countS +=1;
                break;
                }
        }
        if (buf[0] == 'I'){
            for(int i=1; i<150; i++){
                if((buf[i] =='\n') || (buf[i] == ' ')) break;
                tmpbuf[i-1] = buf[i];
            }
            nobs = atoi(tmpbuf);
            memset(tmpbuf,0,sizeof(tmpbuf));
            countI +=1;
            
        }
    }

    if (feof(fp))
       fclose(fp);
}