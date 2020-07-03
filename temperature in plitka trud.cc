#include <stdio.h>
#include <iostream>
#include <stdlib.h>
using namespace std;

#define X 7
#define Y 4
#define T 25

#define dx 0.25
#define dy 0.25
#define dt 0.5
#define k 0.1

int main() {
	int i, j, t;
	int m; double d;
    int xp = X/dx+1, yp=Y/dy+1, tp=T/dt;
    
    double *XX = new double[xp*yp];
    double **Z = new double*[yp];
    for (int i = 0; i < yp; ++i) {
        Z[i]=new double[xp];
    }
    for (int i = 0; i < yp; i++) {
        for (int j = 0; j < xp; j++) {
            Z[i][j]=0.0;
        }
    }
    double **TT = new double*[xp*yp];
    for (int i = 0; i < xp*yp; ++i) {
        TT[i]=new double[xp*yp+1];
    }
    for (int i = 0; i < yp*xp; i++) {
        for (int j = 0; j < xp*yp+1; j++) {
            TT[i][j]=0.0;
        }
        XX[i]=0;
    }

	FILE *points=fopen("points.txt", "w");
	FILE *script=fopen("scr.gnu", "w");
    for(int i=0; i < xp*yp; i++) {
		if((i==xp+1-yp) ||  (( i%xp == (xp-yp+i/xp) ) && i/xp!=0) ||  ( i%xp == (xp+1-yp+i/xp) && i/xp!=0) && i/xp<(yp-1)) {
			TT[i][i]=1;
			TT[i][xp*yp]=100;
		}
       
        else if(i%xp > (xp+1-yp+i/xp)) {
			TT[i][i]=1;
			TT[i][xp*yp]=-10;
		}
        else if(i>xp*(yp-1) && i!=xp*(yp-2)) { //нижняя
            TT[i][i]=1;
            TT[i][yp*xp]=50;
        }
        
        else if(i%xp==0) { //2 род
            TT[i][i]=1;
            TT[i][i+1]=-1;
            TT[i][xp*yp]=0.*dx/k;
        }
		else if(i<xp && i!=0) { //верхняя стенка
            TT[i][i]=1;
            TT[i][xp*yp]=50;
        }
        
        else {
			TT[i][i]= -(2.*k/(dx*dx) + 2.*k/(dy*dy) + 1./dt);
			TT[i][i-1]=k/(dx*dx);
			TT[i][i+1]=k/(dx*dx);
			TT[i][i-(xp)]=k/(dy*dy);
			TT[i][i+(xp)]=k/(dy*dy);
			TT[i][xp*yp]=-XX[i]/dt;
		}
        
    }
	
	
	for(i=0;i<yp; i++) {
		for(j=0; j<xp; j++) {
			Z[i][j]=XX[i*xp+j];
		}
	}
	
		fprintf(points, "1\n");
		for(i=0; i<yp; i++) {
			for(j=0; j<xp; j++) {
				fprintf(points, "%lf %lf %lf\n", j*dx, i*dy, Z[i][j]);
			}
			fprintf(points, "\n");
		}
		fprintf(points, "\n");
	for(t = 2; t<tp; t++) {
		cout << t << " of " << tp << endl;
		//////Прямой ход////////////////////////////
		for(i = 0 ; i < xp*yp-1; i++){
		  for(j = i+1 ; j < yp*xp ; j++){
			d = (TT[j][i]/TT[i][i]);	
			for(m = i; m < xp*yp+1 ; m++)
					TT[j][m]-=(TT[i][m]*d);
		  } 
		}

		////Обратный ход//////////////////////////////
		XX[xp*yp-1] = TT[yp*xp-1][yp*xp-1] / TT[yp*xp-1][yp*xp];
		XX[yp*xp-1] = 1/XX[yp*xp-1];
		for(i = yp*xp-2  ; i >= 0 ; i--) {
			for(j = 0 , d = 0 ; j < xp*yp-i-1 ; j++)
				d+=(XX[yp*xp-j-1]*TT[i][yp*xp-j-1]*(-1));
			XX[i] = (TT[i][yp*xp]+d) / TT[i][i];
		}
		
		for(i=0;i<yp; i++) {
			for(j=0; j<xp; j++) {
				Z[i][j]=XX[i*xp+j];
			}
		}
		
			fprintf(points, "%d\n", t);
			for(i=0; i<yp; i++) {
				for(j=0; j<xp; j++) {
					fprintf(points, "%lf %lf %lf\n", j*dx, (yp-1-i)*dy, Z[i][j]);
				}
				fprintf(points, "\n");
			}
			fprintf(points, "\n");
			for(i=0; i<xp*yp; i++) {
				if(!((i==xp+1-yp) ||  (( i%xp == (xp-yp+i/xp) ) && i/xp!=0) ||  ( i%xp == (xp+1-yp+i/xp) && i/xp!=0) && i/xp<(yp-1))
				&& !(i%xp > (xp+1-yp+i/xp))
				&& !(i>xp*(yp-1) && i!=xp*(yp-2))
				&& !(i%xp==0)
				&& !(i<xp && i!=0))
				TT[i][xp*yp]=-XX[i]/dt;
			}
    }
	fclose(points);
	fprintf(script, "set pm3d map\n");
	//fprintf(script, "set pm3d flush begin ftriangles scansforward interpolate 2,2\n");
	fprintf(script, "set palette\n");
	fprintf(script, "set cbrange [0:200]\n");
	fprintf(script, "set zrange [0:200]\n");
	fprintf(script, "set xrange [0:7]\n");
	fprintf(script, "set yrange [0:4]\n");
	fprintf(script, "do for [i = 0:%d] {\n", tp-2);
	fprintf(script, "\tsplot 'points.txt' index i using 1:2:3 with pm3d\n");
	fprintf(script, "\tpause 0.5\n");
	fprintf(script, "}\npause -1");
	fclose(script);	


	system("gnuplot ./scr.gnu");


}
