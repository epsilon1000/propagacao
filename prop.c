#include <stdio.h>
#include <math.h>
#include <stdlib.h>


//			MODELO
//      _________________________________________  -> (x,i)
//	|
//	|
//	|
//	|				
//	|  |
//	|  V
//	| (z,k)				
//	|	
//	|					0  (y, j)
//	|

//
//			P(i,j,k)=P[Nx*Nz*(j)+Nz*(i)+(k)]
// Es necesario hacer un nuevo for para j (¿?), ademas de emplear fronteras absorbentes para j (y) 
//
//
//
//
#define  A 30		//amplitud
//#define r(i,j,k) 		r[Nx*Nz*(j)+Nz*(i)+(k)]
#define U1(i,j,k)		U1[Nx*Nz*(j)+Nz*(i)+(k)]// Un-1
#define U2(i,j,k)		U2[Nx*Nz*(j)+Nz*(i)+(k)]//Un
#define U3(i,j,k)		U3[Nx*Nz*(j)+Nz*(i)+(k)]//Un+1
#define c(i,j,k)		c[Nx*Nz*(j)+Nz*(i)+k] //velocidad
#define borde 40 

#define dbg(x) fprintf(stderr,"%s \n",x)

float source(float t,float fc);
float laplac3dfd(float *indata, float *outdata, int Nx, int Ny, int Nz);

int main(int argc, char *argv[])
{
int i,j,k,n,ns,nshots,Sx, Profsour;
int nz, Nx, Ny, Nt,amp,salto=10,indice_archivo;
float fc,dx,dy,dz,dt;
char caract[20],velname[20];


ns=0; //************************************************************************************>>>>>>>>>>>> OJO, cambiar para mas disparos
//======================================================================
//			Leyendo Datos del Modelo
//======================================================================

//printf("Leyendo Datos del modelo \n");

//FILE *sintmodel;

Nx=100;
Ny=100;
nz=60;
Nt=1000;
dx=0.01;
dy=0.01;
dz=0.01;
dt=0.002;
fc=8;
nshots=1;
Profsour=41;
/*------------------------------------------------------------------------------------------------------------------------------------------------*/


float *Lsource=calloc(nshots,sizeof(float));

    Lsource[i]=50;  

//======================================================================
//======================================================================
int Nz=nz+borde;    
float factor;
float beta=pow(dt,2)/pow(dx,2);
float *U1 = calloc(Nx*Ny*Nz,sizeof(float));//Un-1
float *U2 = calloc(Nx*Ny*Nz,sizeof(float));//Un
float *outdata = calloc(Nx*Ny*Nz,sizeof(float));//DUn
float *U3 = calloc(Nx*Ny*Nz,sizeof(float));//Un+1
float  *c = calloc(Nx*Ny*Nz,sizeof(float));//velocidad
//float  *r = calloc(Nx*Ny*Nz,sizeof(float));//Ubicacion espacial de la fuente
//======================================================================
//======================================================================
//printf("Leyendo Modelo de Velocidad \n");

//        FILE *modelo;
//        modelo = fopen (velname,"rb");
/***google- >*///	fread(v,Nx*nz*sizeof(float),1,modelo);

for(j=0;j<Ny;j++)
{
	for(i=0;i<Nx;i++)
	{
		for(k=0;k<Nz;k++)
		{
		c(i,j,k)=1;
		}
	}
}	

//fclose(modelo);

//for(i=0;i<Nx;i++){
 //  for(j=0;j<borde;j++){    /*Para que cargar velocidades solo en x?*/
  //     c(i,j)=v(i,0);
 //       c2(i,j)=v(i,0);
 //  }
//}

//for(i=0;i<Nx;i++){
//   for(j=0;j<nz;j++){      /*Al hacer esto no sobreescribimos lo anterior?*/
//       c(i,j+borde)=v(i,j);
//        c2(i,j+borde)=v(i,0);
//   }
       
// }

 
//printf("\n Datos de Modelo de Velocidad cargados...\n");


//======================================================================
//======================================================================
//for(ns=0;ns<nshots;ns++){//loop de los disparos
indice_archivo=0;
//printf("Disparando la fuente %d \n",(ns+1));/*original (ns+1)*/



// for(j=0;j<Ny;j++)
// {
// 	for(i=0;i<Nx;i++)
// 	{
// 		for(k=0;k<Nz;k++)
// 		{
//             printf("%d     %d       %f\n",i,j, outdata[Nx*Nz*(j)+Nz*(i)+(k)] );
// 		}
// 	}
// }	

//====================================================================
//			Programa Principal
//====================================================================
//dbg("Inicia Iteracion en el tiempo");
for(n=0;n<Nt;n++) { //Inicia iteración en el tiempo
//printf("Calculando el campo para el disparo %d tiempo %d \n",ns,n);
//--------------------------------------------------------------------
//Para Obtener Un+1

laplac3dfd(U2,outdata,Nx,Ny,Nz);
// dbg("termina calculo del laplaciano");
for(j=4;j<Ny-4;j++)
{
 for(i=4;i<Nx-4;i++)
 {
   for(k=4;k<Nz-4;k++)
   {
/*=================== Calcula el Un+1 ======================*/
//=========================================================================================================
 
    U3(i,j,k)=2*U2(i,j,k)-U1(i,j,k)+pow(c(i,j,k),2)*beta*outdata[Nx*Nz*(j)+Nz*(i)+(k)];
    U3(50,50,50)+=source(n*dt,fc);	
    if(k==90 && n==500)
		{
			
 		 printf("%d 	%d 	%f \n", i, j, U3(i,j,k));

		}	
   }
 }
}
// dbg("Termina calculo de U3");

//=========================================================================================================
//			Imprimir cada campo en cada tiempo y crear el archivo
//=========================================================================================================

   FILE *archvis;
char  fileoutputname[25];
sprintf(fileoutputname, "forward_%02d_%04d.bin", (ns+1), indice_archivo);
        archvis = fopen(fileoutputname,"wb");
	fwrite(U3,Nx*Ny*Nz*sizeof(float),1,archvis);
fclose(archvis); 
indice_archivo +=1;	    






 //========================================================================================================
//			Actualiza los campos de desplazamiento
//=========================================================================================================
// dbg("Actualizando campos");
for(j=0;j<Ny;j++)
{
 for(i=0;i<Nx;i++) 
  {
   for(k=0;k<Nz;k++)
    {
	U1(i,j,k)=U2(i,j,k);
	U2(i,j,k)=U3(i,j,k);
	U3(i,j,k)=0.0;
    }
  }  
}


} //Finaliza iteracion en el tiempo
/*for(j=0;j<Ny;j++)
{
  for(i=0;i<Nx;i++) 
  {
   for(k=0;k<Nz;k++)
    {
	U1(i,j,k)=0.0;
	U2(i,j,k)=0.0;
	U3(i,j,k)=0.0;
        outdata[Nx*Nz*(j)+Nz*(i)+(k)]=0.0;
       // r(i,j,k)=0.0;
    }
   
  }
}*/

free(U1);free(U2);free(U3);

return 0;
}//END OF MAIN

//========================================================================
float source(float t,float fc)
{
  float R;
  R=A*(1-2*M_PI*M_PI*fc*fc*(t-0.05)*(t-0.05))*exp(-1*M_PI*M_PI*fc*fc*(t-0.05)*(t-0.05));
	return R;
}

//========================================================================
//========================================================================

float laplac3dfd(float *indata,float *outdata, int Nx, int Ny, int Nz)
	
{
  
int i,j,k;

for(j=4;j<Ny-4;j++)
	{	
	for(i=4;i<Nx-4;i++)
		{
		for(k=4;k<Nz-4;k++)
			{
	outdata[Nx*Nz*(j)+Nz*(i)+(k)]=(-1.0/560)*indata[Nx*Nz*(j)+Nz*(i)+(k-4)]
	+(8.0/315)*indata[Nx*Nz*(j)+Nz*(i)+(k-3)]
	+(-1.0/5)*indata[Nx*Nz*(j)+Nz*(i)+(k-2)]+(8.0/5)*indata[Nx*Nz*(j)+Nz*(i)+(k-1)]
	+(-205.0/72)*indata[Nx*Nz*(j)+Nz*(i)+(k)]+(8.0/5)*indata[Nx*Nz*(j)+Nz*(i)+(k+1)]
	+(-1.0/5)*indata[Nx*Nz*(j)+Nz*(i)+(k+2)]+(8.0/315)*indata[Nx*Nz*(j)+Nz*(i)+(k+3)]
	+(-1.0/560)*indata[Nx*Nz*(j)+Nz*(i)+(k+4)];

	outdata[Nx*Nz*(j)+Nz*(i)+(k)]+=(-1.0/560)*indata[Nx*Nz*(j)+Nz*(i-4)+(k)]
	+(8.0/315)*indata[Nx*Nz*(j)+Nz*(i-3)+(k)]
	+(-1.0/5)*indata[Nx*Nz*(j)+Nz*(i-2)+(k)]+(8.0/5)*indata[Nx*Nz*(j)+Nz*(i-1)+(k)]
	+(-205.0/72)*indata[Nx*Nz*(j)+Nz*(i)+(k)]+(8.0/5)*indata[Nx*Nz*(j)+Nz*(i+1)+(k)]
	+(-1.0/5)*indata[Nx*Nz*(j)+Nz*(i+2)+(k)]+(8.0/315)*indata[Nx*Nz*(j)+Nz*(i+3)+(k)]
	+(-1.0/560)*indata[Nx*Nz*(j)+Nz*(i+4)+(k)];

	outdata[Nx*Nz*(j)+Nz*(i)+(k)]+=(-1.0/560)*indata[Nx*Nz*(j-4)+Nz*(i)+(k)]
	+(8.0/315)*indata[Nx*Nz*(j-3)+Nz*(i)+(k)]
	+(-1.0/5)*indata[Nx*Nz*(j-2)+Nz*(i)+(k)]+(8.0/5)*indata[Nx*Nz*(j-1)+Nz*(i)+(k)]
	+(-205.0/72)*indata[Nx*Nz*(j)+Nz*(i)+(k)]+(8.0/5)*indata[Nx*Nz*(j+1)+Nz*(i)+(k)]
	+(-1.0/5)*indata[Nx*Nz*(j+2)+Nz*(i)+(k)]+(8.0/315)*indata[Nx*Nz*(j+3)+Nz*(i)+(k)]
	+(-1.0/560)*indata[Nx*Nz*(j+4)+Nz*(i)+(k)];
	
	//outdata[Nx*Nz*(j)+Nz*(i)+(k)]=2;

			}
		}	
	}	

}
//========================================================================
