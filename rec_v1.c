#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>

typedef uint64_t utype64;
typedef int64_t type64;
typedef __int128 type128;

const utype64 lower58bits = 0x3ffffffffffffff;
const utype64 lower57bits = 0x1ffffffffffffff;


__inline__ utype64 rdtsc() {
   uint32_t lo, hi;
   __asm__ __volatile__ ("xorl %%eax,%%eax \n        cpuid"::: "%rax", "%rbx", "%rcx", "%rdx");
   __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
   return (utype64)hi << 32 | lo;
}

__inline__ utype64 rdtscp() {
   uint32_t lo, hi;
   __asm__ __volatile__ ("rdtscp": "=a" (lo), "=d" (hi) :: "%rcx");
   __asm__ __volatile__ ("cpuid"::: "%rax", "%rbx", "%rcx", "%rdx");
   return (utype64)hi << 32 | lo;
}


/*Toeplitz matrix-vector product where the size is 3 */

void tmvp(type64 *mat, type64 *vec, type128 *result){
	type64 common=mat[0]+mat[1];
	type128 m2=(type128)mat[1]*(vec[0]-vec[1]);
	type128 m3=(type128)mat[0]*(vec[0]-vec[2]);
	type128 m4=(type128)mat[2]*(vec[1]-vec[2]);
	result[0]=m3 + m4 + (type128)(mat[4]+mat[2]+mat[0])*vec[2];
	result[1]=m2 - m4 + (type128)(common+mat[2])*vec[1];
	result[2]=(type128)(mat[3]+common)*vec[0] - m2 - m3;
}

/*Residue multiplication modulo 521-bit Mersenne prime 
FG = Z (mod pow(2,521) - 1) where F and G comprise of 9-limb
F[i] and G[i] are unsigned 58-bit where i=1,2,3,4,5,6,7
F[0], G[0], and Z[0] are at most unsigned 59-bit 
F[8], G[8], and Z[8] are unsigned 57-bit*/

void TMVP_recursive(type64 *F, type64 *G, type64 *Z){
	type64 T1[5]={0,0,0,2*F[6],2*F[5]}, T6[5]; 
	type64 T5=2*F[8], c=2*F[7], temp_vec[3];
	type128 M1[3], M2[3], M3[3], M4[3], M5[3], M6[3];

	T1[0]=G[0]-G[3];		T1[1]=G[1]-G[4];		T1[2]=G[2]-G[5];
	T6[0]=F[3];				T6[1]=F[4];				T6[2]=F[2];				T6[3]=F[5];				T6[4]=F[1];		
	tmvp(T6, T1, M2);

	T1[0]=G[3]-G[6];		T1[1]=G[4]-G[7];		T1[2]=G[5]-G[8];
	T6[0]=T1[3];			T6[1]=c;				T6[2]=T1[4];			T6[3]=T5;				T6[4]=2*F[4];
	tmvp(T6, T1, M4);

	T1[0]=G[0]-G[6];		T1[1]=G[1]-G[7];		T1[2]=G[2]-G[8];
	T6[0]=F[0];				T6[1]=F[1];				T6[2]=T5;				T6[3]=F[2];				T6[4]=c;
	tmvp(T6, T1, M3);

	T6[0]=F[6]+F[3];		T6[1]=F[7]+F[4];		T6[2]=F[5]+F[2];		T6[3]=F[8]+F[5];		T6[4]=F[4]+F[1];	
	T1[0]=T6[0]+F[0];		T1[1]=T6[1]+F[1];		T1[2]=T6[2]+T5;			T1[3]=T6[3]+F[2];		T1[4]=T6[4]+c;
	temp_vec[0]=G[0];		temp_vec[1]=G[1];		temp_vec[2]=G[2];
	tmvp(T1, temp_vec, M1);

	T6[0]+=T1[0];			T6[1]+=T1[1];			T6[2]+=T1[2];			T6[3]+=T1[3];			T6[4]+=T1[4];		
	temp_vec[0]=G[6];		temp_vec[1]=G[7];		temp_vec[2]=G[8];
	tmvp(T6, temp_vec, M6);


	T1[0]+=F[6];			T1[1]=T1[4];			T1[3]=T1[2];			T1[2]=T6[3];			T1[4]=T6[1];
	temp_vec[0]=G[3];		temp_vec[1]=G[4];		temp_vec[2]=G[5];
	tmvp(T1, temp_vec, M5);

	type128 C = M1[2] - M2[2] - M3[2];
	c = ((type64)C)&lower57bits;
	C = (M3[0] + M4[0] + M6[0]) + (C>>57);
	Z[0] = ((type64)C)&lower58bits;
	C = (M3[1] + M4[1] + M6[1]) + (C>>58);
	Z[1] = ((type64)C)&lower58bits;	
	C = (M3[2] + M4[2] + M6[2]) + (C>>58);
	Z[2] = ((type64)C)&lower58bits;
	C = (M2[0] - M4[0] + M5[0]) + (C>>58);
	Z[3] = ((type64)C)&lower58bits;	
	C = (M2[1] - M4[1] + M5[1]) + (C>>58);
	Z[4] = ((type64)C)&lower58bits;	
	C = (M2[2] - M4[2] + M5[2]) + (C>>58);
	Z[5] = ((type64)C)&lower58bits;
	C = (M1[0] - M2[0] - M3[0]) + (C>>58);
	Z[6] = ((type64)C)&lower58bits;
	C = (M1[1] - M2[1] - M3[1]) + (C>>58);
	Z[7] = ((type64)C)&lower58bits; 
	c += (type64)(C>>58);
	Z[8] = c&lower57bits;	
	Z[0] += (c>>57);

}

//////////////////////////////////////////////////////////////////////////////////////////////

int main(){
	FILE *Fptr, *Gptr, *Zptr;
	char  file_F[] = "Fdata.txt", file_G[] = "Gdata.txt", file_Z[] = "Zdata.txt";

	utype64 start, end, min_ccycle=10000;
	type64 F[9], G[9], Z[9];
	unsigned int i, j, jj, k,  values=1000;

	if((Fptr=fopen(file_F, "r"))==NULL || (Gptr=fopen(file_G, "r"))==NULL){
		printf("Cannot open the file(s) for reading.\n");
		return 0;
	}
	if((Zptr=fopen(file_Z, "w"))==NULL ){
			printf("Cannot open the file(s) for writing.\n");
			return 0;
	}
	start=rdtsc();
	end=rdtscp();
	start=rdtsc();
	end=rdtscp();
	for(j=0; j<values; j++){
		for(i=0; i<9; i++){
			fscanf(Fptr,"%lu",&F[i]);
			fscanf(Gptr,"%lu",&G[i]);
		}
		start=rdtsc();
		for (jj=0;jj<1000;jj++)
		{
			TMVP_recursive(F, G, Z);
			TMVP_recursive(Z, G, F);
			//TMVP_recursive(F, Z, G);
			//TMVP_recursive(G, Z, F);
			
		}
		end=rdtscp();
		if((end - start)/2000 < min_ccycle)
			min_ccycle = (end - start)/2000;
		if(j%100 == 0){
			printf("So far, the minimum cycle count is :: %lu\n", min_ccycle);
			for(i=0; i<9; i++)
				fprintf(Zptr, "%lu ", Z[i]);
			fprintf(Zptr, "\n");
		}
	}
	fclose(Fptr);
	fclose(Gptr);
	fclose(Zptr);
	
	printf("The minimum number of clock cycles is :: %lu\n", min_ccycle);
	return 0;
}
