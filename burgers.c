#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

#define PI 4*atan(1)
#define ALPHA 0.5
//#define RANDOM

double* GenerateGaussianWhiteNoise(double sigma, int size);
//Structure of the model
struct ModelBurgers
{
  	double* U_;//burgers variable
	double* Un_;//burgers variable, storage of U previous time
  	double CFL_;//condition CFL	
};

struct ModelBurgersAssimil
{
	double* U_;//burgers variable
	double* Un_;//storage of U, previous time
	double CFL_;//condition CFL
	double lambda_;//assimilation parameter
	int NObs_;//nombre d'observations à faire
	double* ObservationTime_;//array of observation times
	double** ObservationU_;//array of observations
	double** ObservationUNoisy_;//array of observations
};

struct Mesh
{
	int NCells_;/* mesh size */
	double xmin_; /* Xmin */
	double xmax_;/* Xmax*/
	double dx_;/* space step */
	double* Xhalf_;
};

void initMesh(struct Mesh* mesh, int ncells, double xmin, double xmax)
{
	(*mesh).NCells_ = ncells;
	(*mesh).xmin_ = xmin;
	(*mesh).xmax_ = xmax;
	mesh->dx_ = (xmax-xmin)/(ncells);
	mesh->Xhalf_ = (double*)malloc(ncells*sizeof(double));
	int i;
	for(i=0;i<mesh->NCells_;i++)
	{
		mesh->Xhalf_[i] = (i+1.0/2.0)*mesh->dx_;
	}
}

//Init model function
void initModel(struct ModelBurgers *model, int ncells, double* u0)
{
	int i;
	(*model).U_ = (double*)malloc((ncells)*sizeof(double));
	(*model).Un_ = (double*)malloc((ncells)*sizeof(double));
	
	for(i=0;i<ncells;i++)
	{
		(*model).U_[i]=u0[i];
		(*model).Un_[i]=u0[i];
	}
}

//Init model function
void initModelAssimil(struct ModelBurgersAssimil *model, int ncells, int nobs, double* u0, double lambda)
{
	int i,j;
	model->NObs_ = nobs;
	model->U_ = (double*)malloc((ncells)*sizeof(double));
	model->Un_ = (double*)malloc((ncells)*sizeof(double));
	model->ObservationTime_   = (double*) malloc((nobs+1)*sizeof(double) );
	model->ObservationU_      = (double**)malloc((nobs+1)*sizeof(double*));
	model->ObservationUNoisy_ = (double**)malloc((nobs+1)*sizeof(double*));
	model->lambda_ = lambda;
	
	for(i=0;i<nobs+1;i++)
	{
		model->ObservationU_[i] = (double*)malloc(ncells*sizeof(double));
		model->ObservationUNoisy_[i] = (double*)malloc(ncells*sizeof(double));
		model->ObservationTime_[i]=0.0;
	}
	
	for(i=0;i<nobs+1;i++)
	{
		for(j=0;j<ncells;j++)
		{
			model->ObservationU_[i][j] = 0.0;
			model->ObservationUNoisy_[i][j] = 0.0;
		}
	}
	
	for(i=0;i<ncells;i++)
	{
		model->U_[i]=u0[i];
		model->Un_[i]=u0[i];
	}
}
void reinitModelAssimil(struct ModelBurgersAssimil *model, int ncells, double* u0, double lambda)
{
	model->lambda_ = lambda;
	int i;
	
	 for(i=0;i<ncells;i++)
	{
		model->U_[i]=u0[i];
		model->Un_[i]=u0[i];
	}
	 
	//memcpy(model->U_, u0, ncells*sizeof(double));
	//memcpy(model->Un_, u0, ncells*sizeof(double));
}

double ComputeCFL(struct Mesh *mesh)
{
	return mesh->dx_;
}

double ComputeCFLSource(struct ModelBurgersAssimil *model,
						struct Mesh *mesh)
{
	double xi = 0;
	int i;
	for(i=0 ;i < mesh->NCells_; i++)
	{
		xi = MAX(xi,fabs(model->U_[i]));
	}
	return 1.0/(model->lambda_ + xi/mesh->dx_);
}

void ComputeBurgers(struct ModelBurgers *model, struct Mesh *mesh, double dt, int ncells)
{
	//upwind method
	int i;
	for(i=1;i<ncells;i++)
	{
		model->U_[i] = model->Un_[i] - dt/mesh->dx_*(model->Un_[i]*model->Un_[i]/2.0 - model->Un_[i-1]*model->Un_[i-1]/2);
	}
	model->U_[0] = 0;//Dirichlet conditions
	//Godunov method
  
}

void ComputeBurgersAssimil(struct ModelBurgersAssimil *model, struct Mesh *mesh, double dt, int ncells)
{
	//upwind method
	int i;
	/*for(i=1;i<ncells;i++)
	{
		model->U_[i] = model->Un_[i] - dt/mesh->dx_*(model->Un_[i]*model->Un_[i]/2.0 - model->Un_[i-1]*model->Un_[i-1]/2);
	}*/
	//model->U_[0] = 0;//Dirichlet conditions 
	//Engquyst Osher scheme (ce qu'on obtient avec le schéma cinétique)
	for(i=1;i<ncells-1;i++)
	{
		model->U_[i] = model->Un_[i] - dt/mesh->dx_*(  MAX(0,model->Un_[i])*model->Un_[i]/2.0
													 + MIN(0, model->Un_[i+1])*model->Un_[i+1]/2.0
													 - MAX(0,model->Un_[i-1])*model->Un_[i-1]/2
													 - MIN(0, model->Un_[i])*model->Un_[i]/2.0 );
	}
	model->U_[0] = 0;//Dirichlet conditions
	model->U_[ncells-1] = model->U_[ncells-2];//Dirichlet conditions
	
}

void BuildObservations(struct ModelBurgersAssimil *model,
					   int obsnb, double tim, double* u, int ncells)
{
	int i;
	model->ObservationTime_[obsnb-1] = tim;
	for(i=0;i<ncells;i++)
	{
		model->ObservationU_[obsnb-1][i] = u[i];
	}
}

void BuildNoisyObservations(struct ModelBurgersAssimil *model, struct Mesh *mesh,
					   int obsnb, double tim, double* u, int ncells, double sigma, double epsilon)
{
	int i;
	double* noise = NULL;
#ifndef RANDOM
	double x;
	noise = (double*)malloc(ncells*sizeof(double));
	for(i=0;i<ncells;i++)
	{
		x = mesh->xmin_ + ( i+0.5 )*mesh->dx_ ;
		noise[i] = pow(epsilon,(1.0 - ALPHA))*cos(x/epsilon);
		//printf("noise i = %d = %f ,pow = %f\n",i,noise[i],pow(epsilon,(1.0 - ALPHA)));
	}
#else
	noise = GenerateGaussianWhiteNoise(sigma,ncells);
#endif	
	//noise = GenerateGaussianWhiteNoise(sigma,ncells);
	
	model->ObservationTime_[obsnb-1] = tim;
	for(i=0;i<ncells;i++)
	{
		  //printf("noise = %f \n", noise[i]);
		model->ObservationUNoisy_[obsnb-1][i] = u[i]+noise[i];
	}
	free(noise);
}

double* ComputeSourceTermCentral(struct ModelBurgersAssimil *modelassimil,
								double dt, int nobs, int ncells)
{
	double* res;
	int i;
	res=(double*)malloc(ncells*sizeof(double));
	for(i=0;i<ncells;i++)
	{
		res[i] = modelassimil->lambda_*dt*
		(modelassimil->ObservationU_[nobs][i] - modelassimil->Un_[i]);
	}
	return res;
}

double* ComputeSourceTermCentralNoisy(struct ModelBurgersAssimil *modelassimil,
								double dt, int nobs, int ncells)
{
	double* res;
	int i;
	res=(double*)malloc(ncells*sizeof(double));
	for(i=0;i<ncells;i++)
	{
		res[i] = modelassimil->lambda_*dt*
		(modelassimil->ObservationUNoisy_[nobs][i] - modelassimil->Un_[i]);
	}
	return res;
}

double* ComputeSourceTermUpwind(struct ModelBurgersAssimil *modelassimil)
{
	return NULL;
}

void Swap(struct ModelBurgers* model, struct Mesh* mesh)
{
	int i=0;
	for(i=0;i< (*mesh).NCells_;i++)
	{
		model->Un_[i] = model->U_[i];
	}
}

void SwapAssimil(struct ModelBurgersAssimil* model, struct Mesh* mesh)
{
	int i=0;
	for(i=0;i< (*mesh).NCells_;i++)
	{
		model->Un_[i] = model->U_[i];
	}
}

void exportGnuplot(const char *filename, double *xn, double *un, int taille)
{
	FILE* fichier = NULL;
	fichier = fopen(filename, "w+");
	char buffer[512*1000];
	int i;
	int nbchar=0;
	if (fichier != NULL)
	{
		//fprintf(fichier, "#U \t V \t E  \t Domain \n");
		for(i=0;i<taille;i++)
		{
			nbchar+=sprintf(buffer+nbchar, "%f \t %f \n", xn[i], un[i]);
		}
		fprintf(fichier, "%s", buffer);
		fclose(fichier);
	}
	else
	{
		printf("Could not open the file.");
	}
}

void exportObservation(const char* filename, struct ModelBurgersAssimil *model, int ncells, int nobs)
{
	FILE* fichier = NULL;
	fichier = fopen(filename, "w+");
	char buffer[512*1000];
	int i,j;
	int nbchar=0;

	if (fichier != NULL)
	{
		//Première ligne
		for(i=0;i<nobs;i++)
		{
			//fprintf(fichier,"%f \t", model->ObservationTime_[i]);
			nbchar+=sprintf(buffer+nbchar, "%f \t", model->ObservationTime_[i]);
		}
		nbchar+=sprintf(buffer+nbchar, "\n");
		for(j=0;j<ncells;j++)
		{
			for(i=0;i<nobs;i++)
			{
				//fprintf(fichier,"%f \t", model->ObservationU_[i][j])
				nbchar+=sprintf(buffer+nbchar, "%f \t",model->ObservationU_[i][j]);
			}
			nbchar+=sprintf(buffer+nbchar, "\n");
		}

		fprintf(fichier, "%s", buffer);

		fclose(fichier);
	}
	else
	{
		printf("Could not open the file.");
	}

}

void exportNoisyObservation(const char* filename, struct ModelBurgersAssimil *model, int ncells, int nobs)
{
	FILE* fichier = NULL;
	fichier = fopen(filename, "w+");
	char buffer[512*1000];
	int i,j;
	int nbchar=0;

	if (fichier != NULL)
	{
		//1ere ligne
		for(i=0;i<nobs;i++)
		{
			//fprintf(fichier,"%f \t", model->ObservationTime_[i]);
			nbchar+=sprintf(buffer+nbchar, "%f \t", model->ObservationTime_[i]);
		}
		nbchar+=sprintf(buffer+nbchar, "\n");
		for(j=0;j<ncells;j++)
		{
			for(i=0;i<nobs;i++)
			{
				//fprintf(fichier,"%f \t", model->ObservationU_[i][j])
				nbchar+=sprintf(buffer+nbchar, "%f \t",model->ObservationUNoisy_[i][j]);
			}
			nbchar+=sprintf(buffer+nbchar, "\n");
		}

		fprintf(fichier, "%s", buffer);

		fclose(fichier);
	}
	else
	{
		printf("Could not open the file.");
	}

}
double* GenerateGaussianWhiteNoise(double sigma, int size)
{
	double* a = (double*)malloc(size*sizeof(double));
	double* b = (double*)malloc(size*sizeof(double));
	double* res = (double*)malloc(size*sizeof(double));
	int i;
	//Generate first uniform law
	for(i=0;i<size;i++)
	{
	a[i] = 1.0*rand()/RAND_MAX;
	}
	//Generate second uniform law
	for(i=0;i<size;i++)
	{
	b[i] = 1.0*rand()/RAND_MAX;
	}
	//Generate Gaussian noise, Box Muller method
	for(i=0;i<size;i++)
	{
	res[i] = sigma*sqrt(-2*log(a[i]))*cos(2*PI*b[i]);
	}
	free(a);
	free(b);
  return res;
}


double ComputeL2norm(double* uexact, double* u, double size)
{
  int i;
  double l2norm = 0.0;
  for(i=0;i<size;i++)
  {
    l2norm = l2norm + (uexact[i] - u[i])*(uexact[i] - u[i]);
	  //printf("i = %d and diff =%f , ue = %f et u = %f\n",i, uexact[i] - u[i],uexact[i], u[i]);
  }
  return sqrt(l2norm);
}

double ComputeLinfnorm(double* uexact, double* u, double size)
{
  int i;
  double linfnorm = 0.0;
  for(i=0;i<size;i++)
  {
    linfnorm = MAX(linfnorm, fabs(uexact[i] - u[i]));
  }
  return linfnorm;
}

double ComputeFracSobolevNorm(struct Mesh *mesh, double* uexact, double *u, int size,double gamma)
{
	int i,j;
	double norm = 0.0;
	double x,y;
	for(i=0;i<size;i++)
	{
		x = mesh->xmin_ + (i+0.5)*mesh->dx_;
		for(j=0;j<size;j++)
		{
			if (i!=j)
			{
				y = mesh->xmin_ + (j+0.5)*mesh->dx_;
				//printf("i = %d et j = %d diff1 = %f, diff2 = %f \n",i,j,(uexact[i] - u[i])-(uexact[j] - u[j]),uexact[i] - u[i]);
				norm = norm + 
				((uexact[i] - u[i]) - (uexact[j] - u[j]))*( (uexact[i] - u[i]) - (uexact[j] - u[j]))/
				pow(fabs(x - y),1+2*gamma);
				printf("i = %d, j = %d, norm =%f \n ",i,j,pow(fabs(x - y),1.0+2.0*gamma));
			}
		}
	}
	return sqrt(norm);
}

void exportnorm2lambda(double** lambdaarray, int nlambda, int nexp, const char* filename)
{
	FILE* fichier = NULL;
	fichier = fopen(filename, "w+");
	char buffer[512*10000];
	int i,j;
	int nbchar=0;
	
	if (fichier != NULL)
	{
		for(i=0;i<nexp;i++)
		{
			for(j=0;j<nlambda;j++)
			{
				//fprintf(fichier,"%f \t", model->ObservationU_[i][j])
				nbchar+=sprintf(buffer+nbchar, "%f \t",lambdaarray[i][j]);
			}
			nbchar+=sprintf(buffer+nbchar, "\n");
		}
		
		fprintf(fichier, "%s", buffer);
		
		fclose(fichier);
	}
	else
	{
		printf("Could not open the file.");
	}
}

void exportnorminflambda(double** lambdaarray, int nlambda, int nexp, const char* filename)
{
	FILE* fichier = NULL;
	fichier = fopen(filename, "w+");
	char buffer[512*10000];
	int i,j;
	int nbchar=0;
	
	if (fichier != NULL)
	{
		for(i=0;i<nexp;i++)
		{
			for(j=0;j<nlambda;j++)
			{
				//fprintf(fichier,"%f \t", model->ObservationU_[i][j])
				nbchar+=sprintf(buffer+nbchar, "%f \t",lambdaarray[i][j]);
			}
			nbchar+=sprintf(buffer+nbchar, "\n");
		}
		
		fprintf(fichier, "%s", buffer);
		
		fclose(fichier);
	}
	else
	{
		printf("Could not open the file.");
	}
}

int main()
{
	//Definition des paramètres
	printf("Entering the parameters...");
	const double xmin = 0.0;
	const double xmax = 1.0;
	const int ncells = 100;
	const double eps = 0.0000000001;
	
	const int nbstep = 1000000;
	int nstep = 0;
	const double TFINAL = 2.0;
	const int nobstot = 30;
	double tim = 0.0;
	double dt = 0.0;
	
	double lambda =1;//= 50000.0;
	double nblambda = 1;
	double epsilon = 0.01;
	printf("OK \n");
	
	//Construction des structures
	printf("Building the structures...");
	struct Mesh mymesh;
  	struct ModelBurgers mymodel;
	struct ModelBurgersAssimil mymodelassimil;
	printf("OK \n");
	
	int l;
	int exp;
	int nexperiment = 1;
	double** norm2 = (double**)malloc(nexperiment*sizeof(double*));
	double** norminf = (double**)malloc(nexperiment*sizeof(double*));
	double** normsob = (double**)malloc(nexperiment*sizeof(double*));

	for(exp=0;exp<nexperiment;exp++)
	{
		norm2[exp] = (double*)malloc(nblambda*sizeof(double));
		norminf[exp] = (double*)malloc(nblambda*sizeof(double));
		normsob[exp] = (double*)malloc(nblambda*sizeof(double));

	}
	for(exp=0;exp<nexperiment;exp++)
	{
		for(l=0;l<nblambda;l++)
		{
		norm2[exp][l] = 0.0;
		norminf[exp][l] =0.0;
		normsob[exp][l] =0.0;
		}
	}
	
	//Initialisation des structures
	double* u0 = (double*)malloc(ncells*sizeof(double));
	double* u0assimil = (double*)malloc(ncells*sizeof(double));
	int i;
	for(i = 0; i< ncells;i++)
	{
		u0[i] = 0.0;
	}
	for(i=ncells/8 ; i< ncells/4;i++)
	{
		u0[i] = 1.0;
	}
	
	for(i = 0; i< ncells;i++)
	{
		u0assimil[i] = 0.0;
	}
	for(i=ncells/12 ; i< ncells/6;i++)
	{
		u0assimil[i] = 0.75;
	}
	
	/*for(i=0 ; i<ncells;i++)
	 {
	 float x = i-ncells/6.0;
	 float y = i-ncells/2.0;
	 float sigma = ncells/25.0; 
	 u0[i] = exp(-x*x/(2*sigma*sigma)) + exp(-y*y/(2*sigma*sigma));
	 }*/
	
	int count = 0;
	//Initialisation de la seed pour le bruit blanc
	srand(time(NULL));
	for(exp=0;exp<nexperiment;exp++)//pour chaque expérience
	{
		printf("Experience %d/%d \n", exp+1, nexperiment);
		printf("Initializing the structures...");
		initMesh(&mymesh,ncells,xmin,xmax);
		initModel(&mymodel, ncells, u0);
		initModelAssimil(&mymodelassimil, ncells, nobstot, u0assimil, lambda);
		printf("OK \n");
	
		printf("Time loop Burgers...");
		tim = 0.0;
		dt = ComputeCFL(&mymesh);
		char buff[24];
		int aux=0;
		int aux2=1;
		const char* filename = "burgersdet";
		double obstime = 0.01;
		double obsassimil = (TFINAL)/(mymodelassimil.NObs_+1.0);
		double sigma = 0.01;
		
		while(tim < TFINAL && nstep < nbstep)
		{
			ComputeBurgers(&mymodel, &mymesh, dt, mymesh.NCells_);
			Swap(&mymodel, &mymesh);
			tim = tim + dt;
			nstep++;
			/**Build Observations */
			if(tim > aux2*obsassimil+eps)
			{
				BuildObservations(&mymodelassimil,aux2,tim,mymodel.U_, ncells);
				BuildNoisyObservations(&mymodelassimil,&mymesh,aux2,tim,mymodel.U_, ncells,sigma,epsilon);
				aux2++;
			}
			
			/** Export Gnuplot for several steps (movie) */	  
			if(tim > aux*obstime+eps)
			 {
				aux++;
				sprintf(buff,"%s_%04d.txt",filename, aux);
				exportGnuplot(buff, mymesh.Xhalf_,mymodel.U_,ncells); 
			 }
		}
		mymodelassimil.ObservationTime_[nobstot] = TFINAL;
		//exportObservation("obs.txt", &mymodelassimil, ncells, nobstot+1);
		//exportNoisyObservation("obsnoisy.txt", &mymodelassimil, ncells, nobstot+1);
		printf("OK \n");
	
		printf("Time loop Burgers Assimil...\n");
		count = 0;//compteur de lambda
		for (l=0; l<1; l=l+1)
		{
			lambda = 100;//20 + 1.0*l/2000.0*100.0;
			printf("lambda = %f \n", lambda);
			reinitModelAssimil(&mymodelassimil, ncells, u0assimil, lambda);
			
			aux = 0.0;
			tim = 0.0;
			int nobs = 0;
			double tobs = 0.0;
			dt = ComputeCFLSource(&mymodelassimil, &mymesh);
			
	
			double* source = NULL;
			nstep=0;
			const char* filename2 = "burgersdet2";
			char buff2[24];
			double sortie=0;
			while(tim < TFINAL && nstep < nbstep)
			{
				tobs = mymodelassimil.ObservationTime_[nobs];
				while(tim+dt <= tobs)//tant qu'on n'a pas atteint une nouvelle observation
				{
					ComputeBurgersAssimil(&mymodelassimil, &mymesh, dt, mymesh.NCells_);//on avance burgers normal
					SwapAssimil(&mymodelassimil, &mymesh);
					tim = tim+dt;
					nstep++;
					dt = ComputeCFLSource(&mymodelassimil, &mymesh);
					if(tim > aux*obstime+eps)//on fait une sortie si demandé par obstime
					{
						aux++;
						sprintf(buff2,"%s_%04d.txt",filename2, aux);
						exportGnuplot(buff2, mymesh.Xhalf_,mymodelassimil.U_,ncells); 
					}
					if(tim > 0.8*TFINAL && sortie==0)
					{
						norm2[exp][count] = ComputeL2norm(mymodel.U_, mymodelassimil.U_, ncells);
						//printf("La norme L2 au temps %f est %f \n", tim,norm2[exp][count]);
			
						norminf[exp][count] = ComputeLinfnorm(mymodel.U_, mymodelassimil.U_, ncells);
						//printf("La norme inf au temps %f est %f \n", tim,norminf[exp][count]);
						
						normsob[exp][count] = ComputeFracSobolevNorm(&mymesh, mymodel.U_, mymodelassimil.U_, ncells,0.125);

						sortie=1;
					}
				}
				if(nobs < mymodelassimil.NObs_)
				{
					//on avance jusqu'à l'observation
					dt = tobs - tim;
					ComputeBurgersAssimil(&mymodelassimil, &mymesh, dt, mymesh.NCells_);
					SwapAssimil(&mymodelassimil, &mymesh);
					nstep++;
					tim=tobs;
					//on compute l'observation
					dt = ComputeCFLSource(&mymodelassimil, &mymesh);
					dt = MIN(dt, mymodelassimil.ObservationTime_[nobs+1]-tobs);//on ne va pas plus loin que la prochaine observation
					source = ComputeSourceTermCentralNoisy(&mymodelassimil,
										dt,nobs,ncells);
					ComputeBurgersAssimil(&mymodelassimil, &mymesh, dt, mymesh.NCells_);
					for(i=0;i<mymesh.NCells_;i++)
					{
						mymodelassimil.U_[i]=mymodelassimil.U_[i]+source[i];
					}
					SwapAssimil(&mymodelassimil, &mymesh);
					tim = tim+dt;
					nstep++;
					nobs++;
					if(tim > aux*obstime+eps)
					{
						aux++;
						sprintf(buff2,"%s_%04d.txt",filename2, aux);
						exportGnuplot(buff2, mymesh.Xhalf_,mymodelassimil.U_,ncells); 
					}
					free(source);
				}
				else 
				{
						break;
				}
			}
			count =count+1;
		}//fin de la boucle sur les lambda
		free(mymodelassimil.U_);
		free(mymodelassimil.Un_);
		free(mymodel.U_);
		free(mymodel.Un_);
		free(mymodelassimil.ObservationTime_);
		for(i=0;i<nobstot+1;i++)
		{
			free(mymodelassimil.ObservationU_[i]);
			free(mymodelassimil.ObservationUNoisy_[i]);
		}
		free(mymodelassimil.ObservationU_);
		free(mymodelassimil.ObservationUNoisy_);
	}//fin de la boucle sur les expériences
	free(u0);
	free(u0assimil);
	
	exportnorm2lambda(norm2, nblambda, nexperiment,"norm2.txt");
	exportnorminflambda(norminf, nblambda, nexperiment, "norminf.txt");
	exportnorm2lambda(normsob, nblambda, nexperiment,"normsobesp9.txt");

	for(i=0;i<nexperiment;i++)
	{
		free(norm2[i]);
		free(norminf[i]);
		free(normsob[i]);

	}
	free(norm2);
	free(norminf);
	free(normsob);
	free(mymesh.Xhalf_);
	
  	return EXIT_SUCCESS;
}
