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

/*Residue multiplication modulo 521-bit Mersenne prime 
FG = Z (mod pow(2,521) - 1) where F and G comprise of 9-limb
F[i] and G[i] are unsigned 58-bit where i=1,2,3,4,5,6,7
F[0], G[0], and Z[0] are at most unsigned 59-bit 
F[8], G[8], and Z[8] are unsigned 57-bit*/

void TMVP_mixed(type64 *F, type64 *G, type64 *Z){
	type64 T1[5]={0,0,0,2*F[6],2*F[5]}, T6[5]; 
	type64 T5=2*F[8], c=2*F[7];

	T1[0]=G[0]-G[3];		T1[1]=G[1]-G[4];		T1[2]=G[2]-G[5];
	type64 intra_common=F[4]+F[3];
	type64 inter_common=F[1]+F[2];
	type128 m2 = (type128)F[4]*(T1[0]-T1[1]);
    type128 m3 = (type128)F[3]*(T1[0]-T1[2]);
    type128 m4 = (type128)F[2]*(T1[1]-T1[2]);
    type128 X0 = m3 + m4 + (type128)(inter_common+F[3])*T1[2];
    type128 X1 = m2 - m4 + (type128)(intra_common+F[2])*T1[1];
    type128 X2 = (type128)(F[5]+intra_common)*T1[0] - m2 - m3;
    
    T1[0]=G[0]-G[6];		T1[1]=G[1]-G[7];		T1[2]=G[2]-G[8];
    intra_common=T5+F[0];
    m2 = (type128)F[1]*(T1[0]-T1[1]);
    m3 = (type128)F[0]*(T1[0]-T1[2]);
    m4 = (type128)T5*(T1[1]-T1[2]);
    type128 X3 = m3 + m4 + (type128)(intra_common+c)*T1[2];
    type128 X4 = m2 - m4 + (type128)(intra_common+F[1])*T1[1];
    type128 X5 = (type128)(F[0]+inter_common)*T1[0] - m2 - m3;
    
    T1[0]=G[3]-G[6];		T1[1]=G[4]-G[7];		T1[2]=G[5]-G[8];
    intra_common=T1[3]+c;
    m2 = (type128)c*(T1[0]-T1[1]);
    m3 = (type128)T1[3]*(T1[0]-T1[2]);
    m4 = (type128)T1[4]*(T1[1]-T1[2]);
    type128 X6 = m3 + m4 + (type128)(T1[4]+T1[3]+(2*F[4]))*T1[2];
    type128 X7 = m2 - m4 + (type128)(intra_common+T1[4])*T1[1];
    type128 X8 = (type128)(T5+intra_common)*T1[0] - m2 - m3;
    
    T6[0]=F[4]+F[1];		T6[1]=F[5]+F[2];		T6[2]=F[6]+F[3];
    T6[3]=F[7]+F[4];		T6[4]=F[8]+F[5];	
	T1[0]=T6[0]+c;			T1[1]=T6[1]+T5;			T1[2]=T6[2]+F[0];
	T1[3]=T6[3]+F[1];		T1[4]=T6[4]+F[2];
	T6[0]+=T1[0];			T6[1]+=T1[1];			T6[2]+=T1[2];
	T6[3]+=T1[3];			T6[4]+=T1[4];		
	T5=T1[2]+F[6];

	type128 C = ((type128)T1[2]*G[2])+((type128)T1[3]*G[1])+((type128)T1[4]*G[0]) - X2 - X5;
	c = ((type64)C)&lower57bits;
	C = (((type128)T6[0]*G[8])+((type128)T6[1]*G[7])+((type128)T6[2]*G[6]) + X3 + X6) + (C>>57);
	Z[0] = ((type64)C)&lower58bits;
	C = (((type128)T6[1]*G[8])+((type128)T6[2]*G[7])+((type128)T6[3]*G[6]) + X4 + X7) + (C>>58);
	Z[1] = ((type64)C)&lower58bits;	
	C = (((type128)T6[2]*G[8])+((type128)T6[3]*G[7])+((type128)T6[4]*G[6]) + X5 + X8) + (C>>58);
	Z[2] = ((type64)C)&lower58bits;
	C = (((type128)T6[3]*G[5])+((type128)T6[4]*G[4])+((type128)T5*G[3]) + X0 - X6) + (C>>58);
	Z[3] = ((type64)C)&lower58bits;	
	C = (((type128)T6[4]*G[5])+((type128)T5*G[4])+((type128)T1[0]*G[3]) + X1 - X7) + (C>>58);
	Z[4] = ((type64)C)&lower58bits;	
	C = (((type128)T5*G[5])+((type128)T1[0]*G[4])+((type128)T1[1]*G[3]) + X2 - X8) + (C>>58);
	Z[5] = ((type64)C)&lower58bits;
	C = (((type128)T1[0]*G[2])+((type128)T1[1]*G[1])+((type128)T1[2]*G[0]) - X0 - X3) + (C>>58);
	Z[6] = ((type64)C)&lower58bits;
	C = (((type128)T1[1]*G[2])+((type128)T1[2]*G[1])+((type128)T1[3]*G[0]) - X1 - X4) + (C>>58);
	Z[7] = ((type64)C)&lower58bits; 
	c += ((type64)(C>>58));
	Z[8] = c&lower57bits;	
	Z[0] += (c>>57);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

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
			TMVP_mixed(F, G, Z);
			TMVP_mixed(Z, G, F);
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
