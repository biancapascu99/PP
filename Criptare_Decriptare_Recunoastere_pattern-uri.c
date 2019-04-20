#include <stdio.h>
#include <inttypes.h>
#include<stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    unsigned int inaltime;
    unsigned int latime;
} dim;

typedef struct{
    int x, y, inaltime, latime, cifra;
    double cor;
} detectie;

dim dimensiune_imagine(char* nume_fisier)
{
    FILE* f;
    dim d;
   unsigned int dim_img, latime_img, inaltime_img;

   f = fopen(nume_fisier, "rb");
   if(f == NULL)
   	{
   		printf("nu am gasit imaginea sursa din care citesc");
   		return;
   	}

   fseek(f, 18, SEEK_SET);
   fread(&latime_img, sizeof(unsigned int), 1, f);
   fread(&inaltime_img, sizeof(unsigned int), 1, f);
   fclose(f);

   d.inaltime = inaltime_img;
   d.latime = latime_img;
   return d;
}

typedef struct
{
    unsigned char B, G, R;
}   pixel;

void afisare_pixel ( pixel p )
{
    printf("(%u, %u, %u) ", p.R, p.G, p.B );
}

uint32_t xorshift32(uint32_t state[static 1])
{
	uint32_t x = state[0];
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	state[0] = x;
	return x;
}

uint32_t* generare_xorshift32( int R, int start )
{
    uint32_t a[1], * v;
    int i;

    a[0] = start;
    v = (uint32_t*)malloc(R * sizeof(uint32_t));
    for ( i = 1; i < R; i++ )
    {
        a[0] = xorshift32(a);
        v[i] = a[0];
    }
    return v;
}
void permuta_pixeli ( pixel* v, int* p, int dim )
{
    pixel* aux;
    int i, j;
    aux = (pixel*)malloc( dim * sizeof(pixel));
    for ( i = 0; i < dim; i++ )
        aux[i] = v[i];

    for ( i = 0; i < dim; i++ )
        v[p[i]] = aux[i];
}

pixel pixel_xor_uint(pixel p, unsigned int x)
{
    unsigned char* cc = &x;
    pixel p_aux;

    p_aux.B = p.B ^ (*cc);
    p_aux.G = p.G ^ (*(cc + 1));
    p_aux.R = p.R ^ (*(cc + 2));

    return p_aux;
}

pixel pixel_xor_pixel(pixel p1, pixel p2)
{
    pixel p_aux;

    p_aux.B = p1.B ^ p2.B;
    p_aux.G = p1.G ^ p2.G;
    p_aux.R = p1.R ^ p2.R;

    return p_aux;
}


void Durstenfeld_Xorshift(int* p, int n, uint32_t* R)
{
    unsigned int r, k, i;
    int aux;

    for ( k = n-1; k >= 1; k-- )
    {
        r = R[n-k]%(k+1);
        aux = p[r];
        p[r] = p[k];
        p[k] = aux;

    }
    return p;
}

void permutare_inversa(int* p, int n)
{
    int* aux, i;
    aux = (int*)malloc(n * sizeof(int));
    for ( i = 0; i < n; i++ )
        aux[i] = p[i];
    for ( i = 0; i < n; i++ )
        p[aux[i]] = i;
}

pixel* liniarizare_imagine(char* nume_fisier)
{
    FILE* f;
    unsigned int dim_img, latime_img, inaltime_img;
    unsigned char pRGB[3], header[54], aux;
    int padding;

    f = fopen(nume_fisier, "rb+");
    if (f == NULL)
   	{
   		printf("nu am gasit imaginea sursa din care citesc");
   		return NULL;
   	}

    fseek(f, 18, SEEK_SET);
   	fread(&latime_img, sizeof(unsigned int), 1, f);
    fread(&inaltime_img, sizeof(unsigned int), 1, f);
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    fseek(f, 54, SEEK_SET);
	int i,j;
	pixel** P;
	P = (pixel**)malloc(inaltime_img * sizeof(pixel*));
	for(i = 0; i < inaltime_img; i++)
	{
	    P[i] = (pixel*)malloc(latime_img * sizeof(pixel));
		for(j = 0; j < latime_img; j++)
		{
			fread(pRGB, 3, 1, f);

			P[i][j].R = pRGB[2];
			P[i][j].G = pRGB[1];
			P[i][j].B = pRGB[0];
		}
		fseek(f,padding,SEEK_CUR);
	}
	fclose(f);
//printf("\n");
    pixel* M;
    int k = -1;
    M = (pixel*)malloc( inaltime_img * latime_img * sizeof(pixel));
    for ( i = 1; i <= inaltime_img; i++ )
    {
        for ( j = 0; j < latime_img; j++ )
        {
            k++;
            M[k] = P[inaltime_img - i][j];
        }
    }

    for(i = 0; i < inaltime_img; i++)
	    free(P[i]);
    free(P);

    return M;
}


void incarca_imagine_liniarizata( char* nume_imagine )
{
    FILE* fin, * fout;
    dim d;
    int dimensiune, i;
    d = dimensiune_imagine(nume_imagine);
    dimensiune = d.inaltime * d.latime;
    pixel* M;
    M = (pixel*)malloc( dimensiune * sizeof(pixel));
    M = liniarizare_imagine(nume_imagine);

	unsigned char c;
    fin = fopen(nume_imagine, "rb");
    char s[] = "liniarizare_";
    char nume_liniarizat[25];
    strcpy(nume_liniarizat, s);
    strcat(nume_liniarizat, nume_imagine);
    fout = fopen(nume_liniarizat, "wb+");

    fseek(fout, 0, SEEK_SET);

	while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}
    int padding;
    if(d.latime % 4 != 0)
        padding = 4 - (3 * d.latime) % 4;
    else
        padding = 0;

    fseek(fout, 54, SEEK_SET);

    for ( i = 0; i < dimensiune; i++ )
    {
        fwrite(&M[i], sizeof(pixel), 1, fout);
        if ( i % d.latime == 0 && i != 0 )
            fseek(fout,padding,SEEK_CUR);
        fflush(fout);
    }

    fclose(fin);
    fclose(fout);
    free(M);
}

pixel* criptare( char* nume_sursa, char* nume_destinatie, char* fisier_cheie )
{
    int i, dimensiune, * perm;
    unsigned int cheie1, cheie2;
    unsigned int* R;
    pixel* M, * C;
    dim d;
    FILE* fin, * fout, * ftext;

    ftext = fopen(fisier_cheie, "r");
    fscanf(ftext, "%u%u", &cheie1, &cheie2);

	M = liniarizare_imagine(nume_sursa);
	d = dimensiune_imagine(nume_sursa);

	dimensiune = d.inaltime * d.latime;

    R = (unsigned int*)malloc(2 * dimensiune * sizeof(unsigned int));
    R = generare_xorshift32( 2 * dimensiune, cheie1 );

    perm = (int*)malloc(dimensiune * sizeof(int));
    for ( i = 0; i < dimensiune; i++ )
        perm[i] = i;

    Durstenfeld_Xorshift(perm, dimensiune, R);
    permuta_pixeli( M, perm, dimensiune );

    C = (pixel*)malloc(dimensiune * sizeof(pixel));
    C[0] = pixel_xor_uint( pixel_xor_uint( M[0], cheie2 ), R[dimensiune] );
    for ( i = 1; i < dimensiune; i++ )
        C[i] = pixel_xor_uint( pixel_xor_pixel( C[i-1], M[i] ), R[dimensiune + i] );

	unsigned char c;
    fin = fopen(nume_sursa, "rb");
    fout = fopen(nume_destinatie, "wb+");

    fseek(fout, 0, SEEK_SET);

	while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}

	int padding;
    if(d.latime % 4 != 0)
        padding = 4 - (3 * d.latime) % 4;
    else
        padding = 0;

    fseek(fout, 54, SEEK_SET);

    int k = -1, j;
    pixel* CC;
    CC = (pixel*)malloc(dimensiune * sizeof(pixel));

    for ( i = 0; i < dimensiune; i = i + d.latime )
        for ( j = i + d.latime - 1; j >= i; j-- )
        {
            k++;
            CC[k] = C[j];
        }

    for ( i = dimensiune - 1; i >= 0; i-- )
    {
        fwrite(&CC[i], sizeof(pixel), 1, fout);
        if ( i % d.latime == 0 && i != 0 )
            fseek(fout,padding,SEEK_CUR);
        fflush(fout);
    }

    fclose(fin);
    fclose(fout);
    fclose(ftext);
    free(R);
    free(M);
    free(perm);
    return C;
}

pixel* decriptare( char* nume_sursa, char* nume_destinatie, char* fisier_cheie )
{
    int i, cheie1, cheie2, dimensiune, * perm;
    unsigned int* R;
    pixel* M, * C_prim;
    dim d;
    FILE* fin, * fout, * ftext;

    ftext = fopen(fisier_cheie, "r");
    fscanf(ftext, "%d%d", &cheie1, &cheie2);

	M = liniarizare_imagine(nume_sursa);
	d = dimensiune_imagine(nume_sursa);

	dimensiune = d.inaltime * d.latime;

    R = (unsigned int*)malloc(2 * dimensiune * sizeof(unsigned int));
    R = generare_xorshift32( 2 * dimensiune, cheie1 );

    perm = (int*)malloc(dimensiune * sizeof(int));
    for ( i = 0; i < dimensiune; i++ )
        perm[i] = i;

    Durstenfeld_Xorshift(perm, dimensiune, R);
    permutare_inversa(perm, dimensiune);

    C_prim = (pixel*)malloc( dimensiune * sizeof(pixel) );
    C_prim[0] = pixel_xor_uint( pixel_xor_uint( M[0], cheie2 ), R[dimensiune] );
    for ( i = 1; i < dimensiune; i++ )
        C_prim[i] = pixel_xor_uint( pixel_xor_pixel( M[i-1], M[i] ), R[dimensiune + i] );

    permuta_pixeli(C_prim, perm, dimensiune);


    int k = -1, j;
    pixel* CC;
    CC = (pixel*)malloc(dimensiune * sizeof(pixel));

    for ( i = 0; i < dimensiune; i = i + d.latime )
        for ( j = i + d.latime - 1; j >= i; j-- )
        {
            k++;
            CC[k] = C_prim[j];
        }

    fin = fopen(nume_sursa, "rb");
    fout = fopen(nume_destinatie, "wb+");
    fseek(fin,0,SEEK_SET);
    fseek(fout, 0, SEEK_SET);
    unsigned char c;
	while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}

    int padding;
    if(d.latime % 4 != 0)
        padding = 4 - (3 * d.latime) % 4;
    else
        padding = 0;

    fseek(fout, 54, SEEK_SET);

	for ( i = dimensiune - 1; i >= 0; i-- )
    {
        fwrite(&CC[i], sizeof(pixel), 1, fout);
        if ( i % d.latime == 0 && i != 0 )
            fseek(fout,padding,SEEK_CUR);
        fflush(fout);
    }


    free(R);
    free(M);
    free(perm);
    fclose(fin);
    fclose(fout);
    fclose(ftext);
    return C_prim;
}

void test_chi_patrat(char* nume_fisier)
{
    int dim_img, inaltime_img, latime_img, i, j,
        aparitiiR[256], aparitiiG[256], aparitiiB[256];
    double medie, sumaR, sumaG, sumaB;
    unsigned char pRGB[3];
    pixel** P;
    FILE* f;

    f = fopen(nume_fisier, "rb");
    if (f == NULL)
   	{
   		printf("nu am gasit imaginea sursa din care citesc");
   		return;
   	}

    fseek(f, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, f);
    fread(&inaltime_img, sizeof(unsigned int), 1, f);

   	int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

	fseek(f, 54, SEEK_SET);
	for ( i = 0; i < 256; i++ )
	{
	    aparitiiR[i] = 0;
	    aparitiiG[i] = 0;
	    aparitiiB[i] = 0;
	}

	P = (pixel**)malloc(inaltime_img * sizeof(pixel*));
	for(i = 0; i < inaltime_img; i++)
	{
	    P[i] = (pixel*)malloc(latime_img * sizeof(pixel));
		for(j = 0; j < latime_img; j++)
		{
			fread(pRGB, 3, 1, f);

			aparitiiR[pRGB[2]]++;
			aparitiiG[pRGB[1]]++;
			aparitiiB[pRGB[0]]++;
		}
		fseek( f, padding, SEEK_CUR );
	}

    medie = (double)(inaltime_img * latime_img)/(double)256;

    sumaR = 0;
    sumaG = 0;
    sumaB = 0;

    for ( i = 0; i < 256; i++ )
    {
        sumaR += ( ((double)aparitiiR[i] - medie) * ((double)aparitiiR[i] - medie) )/medie;
        sumaG += ( ((double)aparitiiG[i] - medie) * ((double)aparitiiG[i] - medie) )/medie;
        sumaB += ( ((double)aparitiiB[i] - medie) * ((double)aparitiiB[i] - medie) )/medie;
    }
    printf ("(%.2f, %.2f, %.2f)\n", sumaR, sumaG, sumaB);
    fclose(f);
}

void grayscale_image(char* nume_fisier_sursa, char* nume_fisier_destinatie)
{
   FILE *fin, *fout;
   unsigned int dim_img, latime_img, inaltime_img;
   unsigned char pRGB[3], header[54], aux;

   fin = fopen(nume_fisier_sursa, "rb");
   if(fin == NULL)
   	{
   		printf("nu am gasit imaginea sursa din care citesc");
   		return;
   	}

   fout = fopen(nume_fisier_destinatie, "wb+");

   fseek(fin, 2, SEEK_SET);
   fread(&dim_img, sizeof(unsigned int), 1, fin);

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img, sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);

	fseek(fin,0,SEEK_SET);
	unsigned char c;
	while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}
	fclose(fin);

	int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

	fseek(fout, 54, SEEK_SET);
	int i,j;
	for(i = 0; i < inaltime_img; i++)
	{
		for(j = 0; j < latime_img; j++)
		{
			fread(pRGB, 3, 1, fout);
			aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
			pRGB[0] = pRGB[1] = pRGB[2] = aux;
        	fseek(fout, -3, SEEK_CUR);
        	fwrite(pRGB, 3, 1, fout);
        	fflush(fout);
		}
		fseek(fout,padding,SEEK_CUR);
	}
	fclose(fout);
}

unsigned char** extragere_matrice (char* nume_imagine)
{
    FILE* fin;
    int latime_img, inaltime_img;
    unsigned char pRGB[3];
    fin = fopen(nume_imagine, "rb");
    if(fin == NULL)
   	{
   		printf("nu am gasit imaginea sursa din care citesc");
   		return;
   	}
    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);

    int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

	fseek(fin, 54, SEEK_SET);
	int i,j;
	unsigned char** P;
	P = (unsigned char**)malloc(inaltime_img * sizeof(unsigned char*));
	for(i = 0; i < inaltime_img; i++)
	{
	    P[i] = (unsigned char*)malloc(latime_img * sizeof(unsigned char));
		for(j = 0; j < latime_img; j++)
		{
			fread(pRGB, 3, 1, fin);
			P[i][j] = pRGB[0];
		}
		fseek(fin,padding,SEEK_CUR);
	}
	fclose(fin);
	return P;
}

float corr ( unsigned char** S, unsigned char** I, int x, int y,
           int inaltime_sab, int latime_sab, int inaltime_imag, int latime_imag)
{
    int i, j;
    int s_sum, f_sum;
    float s_bara, s_sum_dev, s_dev,
        f_bara, f_sum_dev, f_dev,
        corelatie = 0.0;
    s_sum = 0;
    f_sum = 0;
    for ( i = 0; i < inaltime_sab; i++ )
    {
        for ( j = 0; j < latime_sab; j++ )
        {
            s_sum = s_sum + S[i][j];
            if ( i + x < inaltime_imag && j + y < latime_imag )
                f_sum = f_sum + I[i+x][j+y];
        }
    }
    s_bara = (float)s_sum/(float)(inaltime_sab * latime_sab);
    f_bara = (float)f_sum/(float)(inaltime_sab * latime_sab);

    s_sum_dev = 0.0;
    f_sum_dev = 0.0;
    for ( i = 0; i < inaltime_sab; i++ )
    {
        for ( j = 0; j < latime_sab; j++ )
        {
            s_sum_dev = s_sum_dev + ( (float)S[i][j] - s_bara ) * ( (float)S[i][j] - s_bara );
            if ( i + x < inaltime_imag && j + y < latime_imag )
                f_sum_dev = f_sum_dev + ( (float)I[i+x][j+y] - f_bara ) * ( (float)I[i+x][j+y] - f_bara );
        }
    }

     s_dev = sqrt( (1.0/(float)(inaltime_sab * latime_sab - 1.0)) * s_sum_dev );
     f_dev = sqrt( (1.0/(float)(inaltime_sab * latime_sab - 1.0)) * f_sum_dev );

     for ( i = 0; i < inaltime_sab; i++ )
        for ( j = 0; j < latime_sab; j++ )
            if ( s_dev * f_dev != 0)
            {
                if ( i + x < inaltime_imag && j + y < latime_imag )
                    corelatie =  corelatie + (1.0/(s_dev * f_dev)) * ((float)I[i+x][j+y] - f_bara) * ((float)S[i][j] - s_bara);
                else
                {
                    corelatie =  corelatie + (1.0/(s_dev * f_dev)) * (0.0 - f_bara) * ((float)S[i][j] - s_bara);
                }
            }

    if ( corelatie != 0 )
        corelatie = (1.0/(float)(inaltime_sab * latime_sab)) * corelatie;
    return corelatie;
}

void conturare( char* nume_fisier, int x, int y, int inalt, int lat, pixel c )
{

    FILE* f = fopen(nume_fisier, "rb+");
    if ( f == NULL )
        printf("nu s-a putut deschide fisierul\n");
    unsigned int latime_img, inaltime_img;
    int padding;

    dim d_imag;
    d_imag = dimensiune_imagine(nume_fisier);
    fseek(f, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, f);
    fread(&inaltime_img, sizeof(unsigned int), 1, f);
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

	fseek(f, 54, SEEK_SET);
    int i, j;
    unsigned char pRGB[2];
	for(i = 0; i < inaltime_img; i++)
	{
		for(j = 0; j < latime_img; j++)
		{// printf("%d %d\n", i, j);
			fread(pRGB, 3, 1, f);
            if ( i == x && j >= y && j <= y + lat )
            {
                pRGB[0] = c.B;
                pRGB[1] = c.G;
                pRGB[2] = c.R;
            }
            if ( i == x + inalt && j >= y && j <= y + lat )
            {
                pRGB[0] = c.B;
                pRGB[1] = c.G;
                pRGB[2] = c.R;
            }
            if ( i >= x && i <= (x + inalt) && j == y )
            {
                pRGB[0] = c.B;
                pRGB[1] = c.G;
                pRGB[2] = c.R;
            }
            if ( i >= x && i <= (x + inalt) && j == y + lat )
            {
                pRGB[0] = c.B;
                pRGB[1] = c.G;
                pRGB[2] = c.R;
            }
        	fseek(f, -3, SEEK_CUR);
        	fwrite(pRGB, 3, 1, f);
        	fflush(f);
		}
		fseek(f,padding,SEEK_CUR);
	}
	fclose(f);
}

int compara_detectii (const void * a, const void * b)
{
    detectie a1, b1;
    a1 = *(detectie*)a;
    b1 = *(detectie*)b;
    if (b1.cor - a1.cor > 0 )
        return 1;
    else
        if (b1.cor - a1.cor == 0 )
            return 0;
    return -1;
}

int se_suprapune( detectie a, detectie b )
{
    if ((a.x + a.inaltime >= b.x && a.x + a.inaltime <= b.x + b.inaltime
                        && a.y + a.latime >= b.y && a.y + a.latime <= b.y + b.latime)

                    ||   (a.y == b.y && a.x + a.inaltime >= b.x && a.x + a.inaltime <= b.x + b.inaltime)
                    ||   (a.x == b.x && a.y + a.latime >= b.y && a.y + a.latime <= b.y + b.latime)
                    ||   (a.y == b.y && b.x + b.inaltime >= a.x && b.x + b.inaltime <= a.x + a.inaltime)
                    ||   (a.x == b.x && b.y + b.latime >= a.y && b.y + b.latime <= a.y + a.latime)

                    ||   (b.x + b.inaltime >= a.x && b.x + b.inaltime <= a.x + a.inaltime
                            && b.y + b.latime >= a.y && b.y + b.latime <= a.y + a.latime)

                    || (a.x + a.inaltime >= b.x && a.x + a.inaltime <= b.x + b.inaltime
                            && b.y + b.latime >= a.y && b.y + b.latime <= a.y + a.latime)

                    || (b.x + b.inaltime >= a.x && b.x + b.inaltime <= a.x + a.inaltime
                            && a.y + a.latime >= b.y && a.y + a.latime <= b.y + b.latime))
                            return 1;
    return 0;

}

void eliminare_non_maxime (detectie* D, int k, double suprapunere_max)
{
    int intersectie, arie1, arie2, reuniune, i, j, L1, L2;
    double suprapunere;
    qsort(D, k, sizeof(detectie), compara_detectii);
    for ( i = 0; i < k-1; i++ )
    {
        for (j = i+1; j < k; j++ )
        {
            if (D[i].x != -1 && D[i].y != -1)
                if ( se_suprapune( D[i], D[j] ) )
                {
                    intersectie = abs(D[i].x + D[i].inaltime - D[j].x) * abs(D[i].y + D[j].latime - D[j].y);
                    arie1 = D[i].inaltime * D[i].latime;
                    arie2 = D[j].inaltime * D[j].latime;
                    reuniune = arie1 + arie2 - intersectie;
                    suprapunere = (double)intersectie/(double)reuniune;
                    if (suprapunere > suprapunere_max)
                        D[j].x = D[j].y = -1;
                }
        }
    }
}

detectie* template_matching( char* nume_imagine, char* nume_sablon, double ps, detectie* D, int* k )
{
    int i, j;
    unsigned char** P_sablon;
    unsigned char** P_img;
    dim d_sablon, d_img;
    float cor;

    P_sablon = extragere_matrice(nume_sablon);
    d_sablon = dimensiune_imagine(nume_sablon);
    P_img = extragere_matrice(nume_imagine);
    d_img = dimensiune_imagine(nume_imagine);

    for ( i = 0; i < d_img.inaltime; i++ )
    {
        for ( j = 0; j < d_img.latime; j++ )
        {
            cor = corr( P_sablon, P_img, i, j, d_sablon.inaltime, d_sablon.latime, d_img.inaltime, d_img.latime );
            if ( cor > ps )
            {
                (*k)++;
                D[*k].x = i;
                D[*k].y = j;
                D[*k].inaltime = d_sablon.inaltime;
                D[*k].latime = d_sablon.latime;
                D[*k].cifra = nume_sablon[5] - 48;
                D[*k].cor = cor;
            }
        }
    }
    if ( *k == -1 )
        *k = 0;
    for ( i = 0; i < d_sablon.inaltime; i ++ )
        free(P_sablon[i]);
    free(P_sablon);
    for ( i = 0; i < d_img.inaltime; i ++ )
        free(P_img[i]);
    free(P_img);
    return D;
}

int main()
{
    // PARTEA 1
    printf("PARTEA 1:\n\n");
    FILE* f1;
    char nume_fisier_sursa[20], nume_fisier_criptat[20],
        nume_fisier_cheie[20], nume_fisier_decriptat[20];

    f1 = fopen("tema_part1.txt", "r");
    fscanf(f1, "%s", nume_fisier_sursa);
    printf("Fisier sursa: %s\n", nume_fisier_sursa);
    fscanf(f1, "%s", nume_fisier_criptat);
    printf("Fisier criptat: %s\n", nume_fisier_criptat);
    fscanf(f1, "%s", nume_fisier_cheie);
    printf("Fisier cheie secreta: %s\n", nume_fisier_cheie);
    fscanf(f1, "%s", nume_fisier_decriptat);
    printf("Fisier decriptat: %s\n", nume_fisier_decriptat);

    incarca_imagine_liniarizata(nume_fisier_sursa);         // incarca imaginea "pe dos"

    printf("\nIncepe criptare...\n");
    criptare(nume_fisier_sursa, nume_fisier_criptat, nume_fisier_cheie);
    printf("Criptare incheiata!\n");
    printf("\nIncepe decriptare...\n");
    decriptare(nume_fisier_criptat, nume_fisier_decriptat, nume_fisier_cheie);
    printf("Decriptare incheiata...\n");

    printf("\nRezultate test chi-patrat:\n");
    test_chi_patrat(nume_fisier_sursa);
    test_chi_patrat(nume_fisier_criptat);
    test_chi_patrat(nume_fisier_decriptat);

    fclose(f1);

    // PARTEA 2
    printf("\n\nPARTEA 2:\n\n");
    FILE* f2;

    int i, j;
    char nume_imagine[20], nume_imagine_grey[20],
        nume_sabloane[10][20], nume_sabloane_grey[10][20];

    f2 = fopen("tema_part2.txt", "r");
    fscanf(f2, "%s", nume_imagine);
    printf("Imagine: %s\n", nume_imagine);
    fscanf(f2, "%s", nume_imagine_grey);
    printf("Imagine gri: %s\n\n", nume_imagine_grey);

    grayscale_image(nume_imagine, nume_imagine_grey);

    for ( i = 0; i < 10; i++ )
    {
        fscanf(f2, "%s", nume_sabloane[i]);
        printf("Sablon %d: %s\n", i, nume_sabloane[i]);
        fscanf(f2, "%s\n", nume_sabloane_grey[i]);
        printf("Sablon %d gri: %s\n\n", i, nume_sabloane_grey[i]);
        grayscale_image(nume_sabloane[i], nume_sabloane_grey[i]);
    }

    detectie* D;
    dim d_img;
    d_img = dimensiune_imagine(nume_imagine);
    D = (detectie*)malloc(d_img.inaltime * d_img.latime * 10 * sizeof(detectie));

    int k = -1;
    printf("Cautare detectii...\n");
    for ( i = 0; i < 10; i++ )
        template_matching(nume_imagine_grey, nume_sabloane_grey[i], 0.5, D, &k);
    printf("Detectii gasite!\n\n");

    D = (detectie*)realloc(D, sizeof(detectie) * k);

    pixel lista_culori[10];
    lista_culori[0] = (pixel){ .R = 255, .G = 0, .B = 0};
    lista_culori[1] = (pixel){ .R = 255, .G = 255, .B = 0};
    lista_culori[2] = (pixel){ .R = 0, .G = 255, .B = 0};
    lista_culori[3] = (pixel){ .R = 0, .G = 255, .B = 255};
    lista_culori[4] = (pixel){ .R = 255, .G = 0, .B = 255};
    lista_culori[5] = (pixel){ .R = 0, .G = 0, .B = 255};
    lista_culori[6] = (pixel){ .R = 192, .G = 192, .B = 192};
    lista_culori[7] = (pixel){ .R = 255, .G = 140, .B = 0};
    lista_culori[8] = (pixel){ .R = 128, .G = 0, .B = 128};
    lista_culori[9] = (pixel){ .R = 128, .G = 0, .B = 0};

    printf("Eliminare non-maxime...\n");
    eliminare_non_maxime(D, k, 0.2);

    dim d_sablon;
    d_sablon = dimensiune_imagine(nume_sabloane[0]);
    printf("Incepere colorare...\n");
    for ( i = 0; i < k; i++ )
        conturare(nume_imagine, D[i].x, D[i].y, d_sablon.inaltime, d_sablon.latime, lista_culori[D[i].cifra] );
    return 0;
}
