#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 10

void yontemSec(int *islem);
void denklemAl(int *der,float katsayilar[]);
void altUstSinirAl(float *alt,float *ust);
void hataAl(float *hata);
float denklemHesapla(int der,float katsayilar[],float x);
void denklemYazdir(int der,float katsayilar[]);
void baslangicDegeriAl(float *x);
void analitikTurevAl(int der,float katsayilar[],float turevKatsayilar[]);
void bisection();
void regulaFalsi();
void newtonRaphson();
void matrisInverse();
void gaussElimination();
void gaussSeidal();
void numericalDerivative();
void simpson();
void trapez();
void gregoryNewton();

int main(){
	
	int islem;
	printf("1- Bisection yontemi\n2- Regula-falsi yontemi\n3- Newton-Raphson yontemi\n4- NxN'lik matrisin tersi\n5- Gauss Eliminasyon yontemi\n6- Gauss Seidal yontemi\n7- Sayisal Turev (merkezi, ileri ve geri fark)\n8- Simpson yontemi\n9- Trapez yontemi\n10- Degisken donusumsuz Gregory Newton Enterpolasyonu");
	yontemSec(&islem);
	
	if(islem == 1){
		bisection();
	}else if(islem == 2){
		regulaFalsi();
	}else if(islem == 3){
		newtonRaphson();
	}else if(islem == 4){
		matrisInverse();
	}else if(islem == 5){
		gaussElimination();
	}else if(islem == 6){
		gaussSeidal();
	}else if(islem == 7){
		numericalDerivative();
	}else if(islem == 8){
		simpson();
	}else if(islem == 9){
		trapez();
	}else if(islem == 10){
		gregoryNewton();
	}
	
	return 0;
}

void bisection(){
	
	int der,iteration=1;
	float alt,ust,orta,hata,katsayilar[MAX],altSonuc,ustSonuc,ortaSonuc;
	
	denklemAl(&der,katsayilar);
	denklemYazdir(der,katsayilar);
	altUstSinirAl(&alt,&ust);
	hataAl(&hata);
	altSonuc = denklemHesapla(der,katsayilar,alt);
	ustSonuc = denklemHesapla(der,katsayilar,ust);
	orta = (alt + ust)/2;
	ortaSonuc = denklemHesapla(der,katsayilar,orta);
	
	if(altSonuc*ustSonuc<0){
		while((ust-alt>=hata)&&(ortaSonuc!=0)){
			printf("\nIterasyon %d: %f",iteration,orta);
			iteration++;
			if(ortaSonuc*ustSonuc<0){
				alt = orta;
				altSonuc = denklemHesapla(der,katsayilar,alt);			
			}else {
				ust = orta;
				ustSonuc = denklemHesapla(der,katsayilar,ust);				
			}
			orta = (ust+alt)/2;
			ortaSonuc = denklemHesapla(der,katsayilar,orta);
		}
		printf("\nKok: %f",orta);
	}else{
		printf("\nBu aralikta kok yoktur veya cift sayida kok vardir!");
	}
}

void regulaFalsi(){
	
	int der,iteration=1;
	float alt,ust,orta,hata,katsayilar[MAX],altSonuc,ustSonuc,ortaSonuc;
	
	denklemAl(&der,katsayilar);
	denklemYazdir(der,katsayilar);
	altUstSinirAl(&alt,&ust);
	hataAl(&hata);
	altSonuc = denklemHesapla(der,katsayilar,alt);
	ustSonuc = denklemHesapla(der,katsayilar,ust);
	orta = ((alt*ustSonuc)-(ust*altSonuc))/(ustSonuc-altSonuc);
	ortaSonuc = denklemHesapla(der,katsayilar,orta);
	
	if(altSonuc*ustSonuc<0){
		while((ust-alt>=hata)&&(ortaSonuc!=0)){
			printf("\nIterasyon %d: %f",iteration,orta);
			iteration++;
			if(ortaSonuc*ustSonuc<0){
				alt = orta;	
				altSonuc = denklemHesapla(der,katsayilar,alt);	
			}else {
				ust = orta;	
				ustSonuc = denklemHesapla(der,katsayilar,ust);		
			}
			orta = ((alt*ustSonuc)-(ust*altSonuc))/(ustSonuc-altSonuc);
			ortaSonuc = denklemHesapla(der,katsayilar,orta);
		}
		printf("\nKok: %f",orta);
	}else{
		printf("\nBu aralikta kok yoktur veya cift sayida kok vardir!");
	}
}

void newtonRaphson(){
	
	int der,iteration=1;
	float alt,ust,hata,katsayilar[MAX],turevKatsayilar[MAX],x1,x2,x1Sonuc,turevSonuc,error,altSonuc,ustSonuc;
	
	denklemAl(&der,katsayilar);
	denklemYazdir(der,katsayilar);
	altUstSinirAl(&alt,&ust);
	hataAl(&hata);
	altSonuc = denklemHesapla(der,katsayilar,alt);
	ustSonuc = denklemHesapla(der,katsayilar,ust);
	
	if(altSonuc*ustSonuc<0){
		do{
			baslangicDegeriAl(&x1);	
		}while(x1<alt||x1>ust);
		
		analitikTurevAl(der,katsayilar,turevKatsayilar);
		printf("\nDenklemin turevi : ");
		denklemYazdir(der-1,turevKatsayilar);
		
		do{
			if(iteration!=1){
				x1 = x2;
			}
			x1Sonuc = denklemHesapla(der,katsayilar,x1);
			turevSonuc = denklemHesapla(der-1,turevKatsayilar,x1);	
			x2 = x1 - (x1Sonuc/turevSonuc);
			printf("\nIterasyon %d: %f",iteration,x2);
			iteration++;
			error = fabs((x2-x1));
		}while((error>=hata) && (x1Sonuc!=0));
		
		printf("\nKok : %f",x2);
			
	}else{
		printf("\nBu aralikta kok yoktur veya cift sayida kok vardir!");
	}
			
}

void matrisInverse(){

	float a[MAX][MAX],oran;
	int i,j,k,n;
	
	printf("Matris boyutunu giriniz: ");
	scanf("%d", &n);
	
	printf("Matris elemanlarini giriniz:\n");
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
		    printf("a[%d][%d] = ",i+1,j+1);
		    scanf("%f", &a[i][j]);
		}
	}
	
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if(i==j){
			   	a[i][j+n] = 1;
		    }else{
			   	a[i][j+n] = 0;
			}
		}
	}
	 
	for(i=0;i<n;i++){
		if(a[i][i] == 0.0){
			printf("Hata!");
			exit(0);
		}
		for(j=0;j<n;j++){
			if(i!=j){
				oran = a[j][i]/a[i][i];
				for(k=0;k<2*n;k++){
				    a[j][k] = a[j][k] - oran*a[i][k];
				}
			}
		}
	}
	 
	for(i=0;i<n;i++){
		for(j=n;j<2*n;j++){
		   	a[i][j] = a[i][j]/a[i][i];
		}
	}
	
	printf("\nMatrisin tersi:\n");
	for(i=0;i<n;i++){
		for(j=n;j<2*n;j++){
		   	printf("%0.3f\t",a[i][j]);
		}
		printf("\n");
	}	
}

void gaussElimination(){
	
	float a[MAX][MAX],x[MAX],oran;
	int i,j,k,n;
	
	printf("Denklem sayisini giriniz: ");
	scanf("%d", &n);
	
	printf("Genisletilmis katsayilar matrisinin elemanlarini giriniz:\n");
	for(i=0;i<n;i++){
		for(j=0;j<n+1;j++){
		    printf("a[%d][%d] = ",i+1,j+1);
		    scanf("%f", &a[i][j]);
		}
	}
	
	for(i=0;i<n-1;i++){
		if(a[i][i] == 0.0){
			printf("Hata!");
			exit(0);
		}
		for(j=i+1;j<n;j++){
			oran = a[j][i]/a[i][i];
			for(k=i;k<n+1;k++){
			  	a[j][k] = a[j][k] - oran*a[i][k];
			}
		}
	}

	x[n-1] = a[n-1][n]/a[n-1][n-1];
	
	for(i=n-2;i>=0;i--){
		x[i] = a[i][n];
		for(j=i+1;j<n;j++){
		  	x[i] = x[i] - a[i][j]*x[j];
		}
		x[i] = x[i]/a[i][i];
	}	
	
	printf("\nUst ucgen matris: \n");
	for(i=0;i<n;i++){
		printf("\n");
		for(j=0;j<n+1;j++){
		    printf("%f\t",a[i][j]);
		}
	}
	
	printf("\n\nKokler:\n");
	for(i=0;i<n;i++){
		printf("x[%d] = %0.3f\n",i, x[i]);
	}
	
}

void gaussSeidal(){
	
	float a[MAX][MAX],x[MAX],oran,tmp,hata;
	int i,j,k,n,iteration=0,flag;
	
	printf("Denklem sayisini giriniz: ");
	scanf("%d", &n);
	
	printf("\nGenisletilmis katsayilar matrisinin elemanlarini giriniz:\n");
	for(i=0;i<n;i++){
		for(j=0;j<n+1;j++){
		    printf("a[%d][%d] = ",i+1,j+1);
		    scanf("%f", &a[i][j]);
		}
	}
	
	printf("\nKoklerin baslangic degerlerini veriniz\n");
	for(i=0;i<n;i++){
		printf("x[%d]: ",i+1);
		scanf("%f", &x[i]);
	}
	
	hataAl(&hata);
	
	for(i=0;i<n;i++){
		float max=a[0][i];
		int maxInd=0;
		for(j=1;j<n;j++){
		    if(a[j][i]>max){
				maxInd=j;
				max=a[j][i];
			}
		}
    	for(k=0;k<n+1;k++){
    		tmp = a[i][k];
    		a[i][k]=a[maxInd][k];
    		a[maxInd][k] = tmp;
		}		
	}
	
	for(i=0;i<n;i++){
		printf("\n");
		for(j=0;j<n+1;j++){
		    printf("%f\t",a[i][j]);
		}
	}		
	printf("\n");
	do{
		iteration++;
		flag=0;
		printf("\nIter. %d:   ",iteration);
		for(i=0;i<n;i++){
			tmp = x[i];
			float sum=0;
			for(j=0;j<n;j++){
				if(i!=j){
					sum+= a[i][j]*x[j];
				}
			}
			x[i]=(a[i][n]-sum)/a[i][i];
			if(fabs(x[i]-tmp)>hata){
				flag=1;
			}
			printf("%f   ",x[i]);	
		}
	}while(flag==1);

}

void numericalDerivative(){
	
	int der;
	float katsayilar[MAX],x,fark,ileriSonuc,geriSonuc,merkeziSonuc;
		
	denklemAl(&der,katsayilar);
	denklemYazdir(der,katsayilar);
	
	printf("\nTurevini almak istediginiz noktayi giriniz: ");
	scanf("%f",&x);
	printf("\nFark degerini giriniz: ");
	scanf("%f",&fark);
	
	merkeziSonuc = ((denklemHesapla(der,katsayilar,x+fark) - denklemHesapla(der,katsayilar,x-fark))/(2*fark));
	ileriSonuc = ((denklemHesapla(der,katsayilar,x+fark) - denklemHesapla(der,katsayilar,x))/fark);
	geriSonuc = ((denklemHesapla(der,katsayilar,x) - denklemHesapla(der,katsayilar,x-fark))/fark);
		
	printf("\nGeri fark turev sonucu: %f",geriSonuc);
	printf("\nMerkezi turev sonucu: %f",merkeziSonuc);
	printf("\nIleri fark turev sonucu: %f",ileriSonuc);
}

void simpson(){
	
	int i,der,n;
	float alt,ust,h,sonuc,katsayilar[MAX],sonuclar[MAX];
	
	denklemAl(&der,katsayilar);
	denklemYazdir(der,katsayilar);
	altUstSinirAl(&alt,&ust);

	do{
		printf("\nAralik sayisini giriniz(cift sayi): ");
		scanf("%d",&n);	
	}while((n % 2 == 1)||n<2);	
	
	h=(ust-alt)/n;
	
	for(i=0;i<n+1;i++){
		sonuclar[i] = denklemHesapla(der,katsayilar,(alt+i*h));
		
		if(i==0){
			sonuc+= sonuclar[i];
		}else if(i==n){
			sonuc+= sonuclar[i];
		}else if(i % 2==1){
			sonuc+= 4*sonuclar[i];
		}else{
			sonuc+= 2*sonuclar[i];
		}
	}
	
	sonuc = (h/3)*sonuc;
	printf("Integral sonucu: %f",sonuc);
	
}

void trapez(){
	
	int i,der,n;
	float alt,ust,h,sonuc,katsayilar[MAX],sonuclar[MAX];

	denklemAl(&der,katsayilar);
	denklemYazdir(der,katsayilar);
	altUstSinirAl(&alt,&ust);
	
	printf("\nAralik sayisini giriniz: ");
	scanf("%d",&n);	
	h=(ust-alt)/n;
	
	for(i=0;i<n+1;i++){
		sonuclar[i] = denklemHesapla(der,katsayilar,(alt+i*h));
		
		if(i==0){
			sonuc+= sonuclar[i];
		}else if(i==n){
			sonuc+= sonuclar[i];
		}else{
			sonuc+= 2*sonuclar[i];
		}
	}
	
	sonuc = (h/2)*sonuc;
	printf("Integral sonucu: %f",sonuc);	
}

void gregoryNewton(){
	
	int i,j=1,k,n,flag=0;
	float A[MAX][MAX],h,hpower=1,fakto=1,x,xSonuc,sonuc;
	
	printf("\nGirilecek deger sayisini giriniz: ");
	scanf("%d",&n);
	printf("\nDegerler arasi farki giriniz: ");
	scanf("%f",&h);
	printf("\nIlk x degerini giriniz: ");
	scanf("%f",&A[0][0]);
	printf("\nf(%.1f) degerini giriniz: ",A[0][0]);
	scanf("%f",&A[0][1]);
	
	for(i=1;i<n;i++){
		printf("\nf(%.1f) degerini girin: ",A[0][0]+i*h);
		scanf("%f",&A[i][1]);
	}
	
	printf("\nHesaplamak istediginiz degeri giriniz: ");
	scanf("%f",&x);
	
	while(flag==0){
		flag=1;
		for(i=0;i<n-j;i++){
			A[i][j+1]= A[i+1][j]-A[i][j];
		}
		
		for(i=0;i<n-j-1;i++){
			if(A[i][j+1]!=A[i+1][j+1]){
				flag=0;
			}
		}
		j++;
	}
	
	xSonuc = A[0][1];
	printf("\nEnterpolasyon denklemi: %.1f",A[0][1]);
	for(i=0;i<j-1;i++){
		
		hpower=1;
		for(k=0;k<i+1;k++){
			hpower = hpower * h;
		}
		
		fakto=1;
		for(k=0;k<i;k++){
			fakto=fakto*(k+2);
		}
		
		sonuc=1;
		printf(" + %.2f",A[0][i+2]/(fakto*hpower));
		for(k=0;k<i+1;k++){
			printf("(x-%.1f)",A[0][0]+(h*k));
			sonuc = sonuc * (x-(A[0][0]+(h*k)));
		}
		xSonuc += sonuc*(A[0][i+2]/(fakto*hpower));
	}
	
	printf("\nf(%.1f): %.2f",x,xSonuc);
}

void yontemSec(int *islem){
	do{
		printf("\n\n1-10 arasinda bir islem seciniz: ");
		scanf("%d",islem);
	}while(*islem<1||*islem>10);
}

void denklemAl(int *der,float katsayilar[]){
	int i;
	printf("\nGirmek istediginiz denklem kacinci dereceden? : ");
	scanf("%d",der);
	
	for(i = *der; i>=0; i--){
		printf("\n%d. dereceden elemanin katsayisini giriniz : ",i);
		scanf("%f",&katsayilar[i]);
	}
	
}

void altUstSinirAl(float *alt,float *ust){
	printf("\nKapali araligin alt sinirini veriniz : ");
	scanf("%f",alt);
	printf("\nKapali araligin ust sinirini veriniz : ");
	scanf("%f",ust);
}

void hataAl(float *hata){
	printf("\nHata farkini veriniz : ");
	scanf("%f",hata);
}

float denklemHesapla(int der,float katsayilar[],float x){
	int i,j;
	float result=0,k;
	for(i = der; i>=0; i--){
		k=1;
		for(j = 0; j<i; j++){
			k=k * x;
		}
		result += k * katsayilar[i];					
	}
	return result;
}

void denklemYazdir(int der,float katsayilar[]){
	int i;
	printf("\n");
	for(i = der; i>=0;i--){
		if(katsayilar[i]!=0){
			printf("%.2fx^%d  ",katsayilar[i],i);
		}
	}
	printf("\n");
}

void baslangicDegeriAl(float *x){
	printf("\nBaslangic degeri veriniz : ");
	scanf("%f",x);
}

void analitikTurevAl(int der,float katsayilar[],float turevKatsayilar[]){
	int i;
	for(i=der;i>0;i--){
		turevKatsayilar[i-1] = katsayilar[i]*i;
	}
}




