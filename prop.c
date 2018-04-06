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
//
//
//
//
#define  A 10		//amplitud
//#define r(i,j,k) 		r[Nx*Nz*(j)+Nz*(i)+(k)]
#define U1(i,j,k)		U1[Nx*Nz*(j)+Nz*(i)+(k)]// Un-1
#define U2(i,j,k)		U2[Nx*Nz*(j)+Nz*(i)+(k)]//Un
#define U3(i,j,k)		U3[Nx*Nz*(j)+Nz*(i)+(k)]//Un+1
#define U12(i,j,k)		U12[Nx*Nz*(j)+Nz*(i)+(k)]//Un+1 2
#define U22(i,j,k)		U22[Nx*Nz*(j)+Nz*(i)+(k)]//Un 2
#define U32(i,j,k)		U32[Nx*Nz*(j)+Nz*(i)+(k)]//Un+1 2
#define c(i,j,k)		c[Nx*Nz*(j)+Nz*(i)+(k)] //velocidad
#define c2(i,j,k)		c2[Nx*Nz*(j)+Nz*(i)+(k)] //velocidad 2
#define seis(i,j,n)		seis[Nx*Nt*(j)+Nt*(i)+(n)] //Sismograma capa 1
#define seis2(i,j,n)		seis2[Nx*Nt*(j)+Nt*(i)+(n)] //Sismograma capa 2
#define seisf(i,j,n)		seis2[Nx*Nt*(j)+Nt*(i)+(n)] //Sismograma diferencia de 1 y 2
#define DDseis(i,n)		DDseis[Nt*(i)+(n)] //Sismograma en 2D para visualizar 
#define DDseis2(i,n)		DDseis2[Nt*(i)+(n)]
#define DDseisf(i,n)		DDseisf[Nt*(i)+(n)]
#define borde	40 
#define alpha	0.008  
#define dbg(x) fprintf(stderr,"%s \n",x)


#define U3_2D(i,j)		U3_2D[Nz*(i)+(k)]//Un+1 2D


float source(float t,float fc);
float laplac3dfd(float *indata, float *outdata, int Nx, int Ny, int Nz);

int main(int argc, char *argv[])
{
int i,j,k,n,ns,nshots,Sx, Profsour;
int nz, Nx, Ny, Nt,amp,salto=10,indice_archivo,indice_archivo2;
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
nz=100;
//Nz=100;
Nt=1000;
dx=0.01;
dy=0.01;
dz=0.01;
dt=0.001;
fc=20;
nshots=1;
Profsour=borde+1;
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
float *c = calloc(Nx*Ny*Nz,sizeof(float));//velocidad
float *seis = calloc(Nx*Ny*Nt,sizeof(float)); //Sismograma capa 1
//------------- Segunda capa -----------------//
float *U12 = calloc(Nx*Ny*Nz,sizeof(float));//Un-1 2
float *U22 = calloc(Nx*Ny*Nz,sizeof(float));//Un 2
float *outdata2 = calloc(Nx*Ny*Nz,sizeof(float));//DUn 2
float *U32 = calloc(Nx*Ny*Nz,sizeof(float));//Un+1 2
float *c2 = calloc(Nx*Ny*Nz,sizeof(float));//velocidad
float *seis2 = calloc(Nx*Ny*Nt,sizeof(float)); //Sismograma capa 2
float *seisf = calloc(Nx*Ny*Nt,sizeof(float)); //Sismograma diferencia 1 y 2
float *DDseis = calloc(Nx*Nt,sizeof(float)); //Sismograma 2D para visualizar 
float *DDseis2 = calloc(Nx*Nt,sizeof(float));
float *DDseisf = calloc(Nx*Nt,sizeof(float));
//float  *r = calloc(Nx*Ny*Nz,sizeof(float));//Ubicacion espacial de la fuente

float	*U3_2D = calloc(Nx*Nz,sizeof(float)); //Un+1 2D
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
		c(i,j,k)=2.0;
		c2(i,j,k)=2.0;
        if(k>=50)
        {
        c(i,j,k)=4.0;
        }
        
		}
	}
}	

// for(j=0;j<Ny;j++)
// {
// 	for(i=0;i<Nx;i++)
// 	{
// 		for(k=10;k<Nz;k++)
// 		{
// 		c(i,j,k)=2;
// 		if(k>50)
// 		   {
// 			c(i,j,k)=4;
// 		   }			   
// 		c2(i,j,k)=1;
// 		}
// 	}
// }	
//fclose(modelo);
 
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
//			¡Programa Principal!
//====================================================================
//dbg("Inicia Iteracion en el tiempo");
//
for(n=0;n<Nt;n++) { //Inicia iteración en el tiempo
	
//printf("Calculando el campo para el disparo %d tiempo %d \n",ns,n);

//--------------------------------------------------------------------






/*=================== Calcula el Un+1 ======================*/
//====================================================================
//Para Obtener Un+1

laplac3dfd(U2,outdata,Nx,Ny,Nz);
// dbg("termina calculo del laplaciano");
for(j=4;j<Ny-4;j++)
{
 for(i=4;i<Nx-4;i++)
 {
   for(k=4;k<Nz-4;k++)
   {
 
    U3(i,j,k)=2*U2(i,j,k)-U1(i,j,k)+pow(c(i,j,k),2)*beta*outdata[Nx*Nz*(j)+Nz*(i)+(k)];
    U3(50,50,Profsour)+=source(n*dt,fc);	//fuente	
    if(j==50 && n==700)

		{
			
// 		 printf("%d 	%d 	%f \n", i, k, U3(i,j,k));

		}	
   }
 }
}
// dbg("Termina calculo de U3");

//=========================================================================================================
//					Fronteras Absorbenters
//=========================================================================================================

//dbg("Comienza el calculo de fronteras absorbentes en x");
//Frontera x=0 ^ x=Nx

    for(j=0;j<Ny;j++)
    {
	for(k=0;k<Nz;k++) 
	{	
	   for(i=0;i<borde;i++) //  0<x<borde
	   {
   		factor = exp(-(alpha*(borde-i))*(alpha*(borde-i)));			   
		U3(i,j,k)=U3(i,j,k)*factor;
		U2(i,j,k)=U2(i,j,k)*factor;
		U1(i,j,k)=U1(i,j,k)*factor;
	   }


	   
	   for(i=Nx-borde;i<Nx;i++) // borde<x<Nx
	   {
   		factor = exp(-(alpha*(Nx-borde-i))*(alpha*(Nx-borde-i)));			   
		U3(i,j,k)=U3(i,j,k)*factor;
		U2(i,j,k)=U2(i,j,k)*factor;
		U1(i,j,k)=U1(i,j,k)*factor;
	   }
	}   
    }	

//dbg("Comienza el calculo de fronteras absorbentes en z");
// Frontera z=0 ^ z=Nz

    for(j=0;j<Ny;j++)
    {
	for(i=0;i<Nx;i++) 
	{	
	  for(k=0;k<borde;k++) //  0<z<borde
	   {
  		factor = exp(-(alpha*(borde-k))*(alpha*(borde-k)));			   
		U3(i,j,k)=U3(i,j,k)*factor;
		U2(i,j,k)=U2(i,j,k)*factor;
		U1(i,j,k)=U1(i,j,k)*factor;
	   }

// Tapa de superior del cubo
	   
	   for(k=Nz-borde;k<Nz;k++) // borde<z<Nz
	   {
   		factor = exp(-(alpha*(Nz-borde-k))*(alpha*(Nz-borde-k)));			   
		U3(i,j,k)=U3(i,j,k)*factor;
		U2(i,j,k)=U2(i,j,k)*factor;
		U1(i,j,k)=U1(i,j,k)*factor;
	   }
	}   
    }	

// //dbg("Comienza el calculo de fronteras absorbentes en y");
//Frontera j=0 ^ j=Ny

    for(i=0;i<Nx;i++)
    {
	for(k=0;k<Nz;k++) 
	{	
	   for(j=0;j<borde;j++) //  0<y<borde
	   {
   		factor = exp(-(alpha*(borde-j))*(alpha*(borde-j)));			   
		U3(i,j,k)=U3(i,j,k)*factor;
		U2(i,j,k)=U2(i,j,k)*factor;
		U1(i,j,k)=U1(i,j,k)*factor;
	   }


	   
	   for(j=Ny-borde;j<Ny;j++) // borde<y<Ny
	   {
   		factor = exp(-(alpha*(Ny-borde-j))*(alpha*(Ny-borde-j)));			   
		U3(i,j,k)=U3(i,j,k)*factor;
		U2(i,j,k)=U2(i,j,k)*factor;
		U1(i,j,k)=U1(i,j,k)*factor;
	   }

	   for(j=0;j<Ny;j++)
           {		   
	
   	        if(j==50)
	        {
                	
		U3_2D(i,k)=U3(i,50,k);
// 		 printf("%d 	%d 	%f \n", i, k, U3_2D(i,k));

		}	
	   }
	}   
    }	



//=========================================================================================================
//			Imprimir cada campo en cada tiempo y crear el archivo
//=========================================================================================================
// //for(j=0;j<Nj;j++) 
//  // {	
        FILE *archvis;
        char  fileout[25];
        sprintf(fileout, "./Datos/forward2D_%02d_%04d.bin", (ns+1), indice_archivo);
        archvis = fopen(fileout,"wb");
        fwrite(U3_2D,Nx*Nz*sizeof(float),1,archvis);
        fclose(archvis); 
        indice_archivo +=1;	    
 
//  // }
// 



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


// //========================================================================================================
// //					Sismograma 
// //========================================================================================================

for(j=0;j<Ny;j++)
{
	for(i=0;i<Nx;i++)
	{
		seis(i,j,n)=U2(i,j,borde);
// //============================================================================================
// //                   Calcula Sismograma en 2D
// //============================================================================================        
		if(j==50)
		{	
			DDseis(i,n)=U2(i,50,borde);
		}


	}	
}

FILE *vel1;
char fileoutputname1[25];
size_t nsize1;
sprintf(fileoutputname1, "./Datos/uno_sismograma_%04d.bin", (ns+1));
vel1 = fopen(fileoutputname1,"w");
nsize1=fwrite(DDseis,Nx*Nt,sizeof(float),vel1);
fclose(vel1);


} //Finaliza iteracion en el tiempo

free(U1);
free(U2);
free(U3);



//==============================================================================================//
//=========== 		Nueva Iteracion para restar los campos (Ui2)		 ===============//
//==============================================================================================//


for(n=0;n<Nt;n++)
{ //Inicia iteración en el tiempo
/*=================== Calcula el Un+1 ======================*/
//====================================================================
//Para Obtener Un+1

	laplac3dfd(U22,outdata2,Nx,Ny,Nz);
// dbg("termina calculo del laplaciano");
	for(j=4;j<Ny-4;j++)
	{
	 for(i=4;i<Nx-4;i++)
 	{
 	  for(k=4;k<Nz-4;k++)
	   {
 
	    U32(i,j,k)=2*U22(i,j,k)-U12(i,j,k)+pow(c2(i,j,k),2)*beta*outdata2[Nx*Nz*(j)+Nz*(i)+(k)];
    	    U32(50,50,Profsour)+=source(n*dt,fc);	
	    if(k==50 && n==600)
		{
			
 //		 printf("%d 	%d 	%f \n", i, j, U32(i,j,k));

		}	
	   }
	 }
       }


//=========================================================================================================
//					Fronteras Absorbenters
//=========================================================================================================

//dbg("Comienza el calculo de fronteras absorbentes en x");
//Frontera x=0 ^ x=Nx

    for(j=0;j<Ny;j++)
    {
	for(k=0;k<Nz;k++) 
	{	
	   for(i=0;i<borde;i++) //  0<x<borde
	   {
   		factor = exp(-(alpha*(borde-i))*(alpha*(borde-i)));			   
		U32(i,j,k)=U32(i,j,k)*factor;
		U22(i,j,k)=U22(i,j,k)*factor;
		U12(i,j,k)=U12(i,j,k)*factor;
	   }


	   
	   for(i=Nx-borde;i<Nx;i++) // borde<x<Nx
	   {
   		factor = exp(-(alpha*(Nx-borde-i))*(alpha*(Nx-borde-i)));			   
		U32(i,j,k)=U32(i,j,k)*factor;
		U22(i,j,k)=U22(i,j,k)*factor;
		U12(i,j,k)=U12(i,j,k)*factor;
	   }
	}   
    }	

//dbg("Comienza el calculo de fronteras absorbentes en z");
// Frontera z=0 ^ z=Nz

    for(j=0;j<Ny;j++)
    {
	for(i=0;i<Nx;i++) 
	{	
	  for(k=0;k<borde;k++) //  0<z<borde
	   {
   		factor = exp(-(alpha*(borde-k))*(alpha*(borde-k)));			   
		U32(i,j,k)=U32(i,j,k)*factor;
		U22(i,j,k)=U22(i,j,k)*factor;
		U12(i,j,k)=U12(i,j,k)*factor;
	   }

// Tapa de superior del cubo
	   
	   for(k=Nz-borde;k<Nz;k++) // borde<z<Nz
	   {
   		factor = exp(-(alpha*(Nz-borde-k))*(alpha*(Nz-borde-k)));			   
		U32(i,j,k)=U32(i,j,k)*factor;
		U22(i,j,k)=U22(i,j,k)*factor;
		U12(i,j,k)=U12(i,j,k)*factor;
	   }
	}   
    }	

//dbg("Comienza el calculo de fronteras absorbentes en y");
//Frontera j=0 ^ j=Ny

    for(i=0;i<Nx;i++)
    {
	for(k=0;k<Nz;k++) 
	{	
	   for(j=0;j<borde;j++) //  0<y<borde
	   {
   		factor = exp(-(alpha*(borde-j))*(alpha*(borde-j)));			   
		U32(i,j,k)=U32(i,j,k)*factor;
		U22(i,j,k)=U22(i,j,k)*factor;
		U12(i,j,k)=U12(i,j,k)*factor;
	   }


	   
	   for(j=Ny-borde;j<Ny;j++) // borde<y<Ny
	   {
   		factor = exp(-(alpha*(Ny-borde-j))*(alpha*(Ny-borde-j)));			   
		U32(i,j,k)=U32(i,j,k)*factor;
		U22(i,j,k)=U22(i,j,k)*factor;
		U12(i,j,k)=U12(i,j,k)*factor;

	   }

	   for(j=0;j<Ny;j++)
          {		   
  	        if(k==50 && n==700) 
	        {

		
//		 printf("%d 	%d 	%f \n", i, j, U32(i,j,k));
            }
            
		}	
	}   
    }
	
       FILE *archvis2;
       char  fileout2[25];
       sprintf(fileout2, "./Datos/DOS_forward2D_%02d_%04d.bin", (ns+1), indice_archivo2);
       archvis2 = fopen(fileout2,"wb");
       fwrite(U3_2D,Nx*Nz*sizeof(float),1,archvis2);
       fclose(archvis2); 
       indice_archivo2 +=1;	    

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
		U12(i,j,k)=U22(i,j,k);
		U22(i,j,k)=U32(i,j,k);
		U32(i,j,k)=0.0;
    	    }
 	  }  
	}


//========================================================================================================
//					Sismograma 2
//========================================================================================================

    for(j=0;j<Ny;j++)
	{	
	
        for(i=0;i<Nx;i++)
		{
			seis2(i,j,n)=U22(i,j,borde);

		if(j==50)
			{	
			DDseis2(i,n)=U22(i,50,borde);
			}
		}
	}		
	

FILE *vel2;
char fileoutputname2[25];
size_t nsize2;
sprintf(fileoutputname2, "./Datos/dos_sismograma_%04d.bin", (ns+1));
vel2 = fopen(fileoutputname2,"w");
nsize2=fwrite(DDseis2,Nx*Nt,sizeof(float),vel2);
fclose(vel2);
//---------------------------------------------------------------------------------------------------------
}//Finaliza iteracion en el tiempo

free(U12);
free(U22);
free(U32);					

//=========================================================================================================
//					RESTA SISMOGRAMAS
//=========================================================================================================


	for(j=0;j<Ny;j++)
	{
		for(i=0;i<Nx;i++)
		{
			
			for(n=0;n<Nt;n++)
			{	
				seisf(i,j,n)=seis(i,j,n)-seis2(i,j,n);
				if(j==50)
				{
				
			//	DDseisf(i,n)=seisf(i,50,borde);
				DDseisf(i,n)=DDseis(i,n)-DDseis2(i,n);
			//	printf("%d	%d	%f \n",i,n,DDseisf(i,n));
				}
			}	
		}	


	}		
   


FILE *vel3;
char fileoutputname3[25];
size_t nsize3;
sprintf(fileoutputname3, "./Datos/fin_sismograma_%04d.bin", (ns+1));
vel3 = fopen(fileoutputname3,"w");
nsize3=fwrite(DDseisf,Nx*Nt,sizeof(float),vel3);
fclose(vel3);
												
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
	

			}
		}	
	}	

}
