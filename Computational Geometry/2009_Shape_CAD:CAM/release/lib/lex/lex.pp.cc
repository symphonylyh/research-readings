# include <cstdio>
extern FILE *fn;
# define U(x) x
# define NLSTATE yyprevious=YYNEWLINE
# define BEGIN yybgin = yysvec + 1 +
# define INITIAL 0
# define YYLERR yysvec
# define YYSTATE (yyestate-yysvec-1)
# define YYOPTIM 1
# define YYLMAX 200
# define output(c) putc(c,yyout)
# define input() (((yytchar=yysptr>yysbuf?U(*--yysptr):getc(fn))==10?(yylineno++,yytchar):yytchar)==EOF?0:yytchar)
# define unput(c) {yytchar= (c);if(yytchar=='\n')yylineno--;*yysptr++=yytchar;}
# define yymore() (yymorfg=1)
# define ECHO fprintf(yyout, "%s",yytext)
# define REJECT { nstr = yyreject(); goto yyfussy;}
int yyleng; extern char yytext[];
int yymorfg;
extern char *yysptr, yysbuf[];
int yytchar;
FILE *yyin = {stdin}, *yyout = {stdout};
extern int yylineno;
struct yysvf { 
	struct yywork *yystoff;
	struct yysvf *yyother;
	int *yystops;};
struct yysvf *yyestate;
extern struct yysvf yysvec[], *yybgin;
extern int yylook();
#ifdef __cplusplus
extern "C" {
#endif
extern int yywrap();
extern int yylex();
extern int yyreject();
extern int yyracc(int);
extern int yyless(int);
#ifdef __cplusplus
}
#endif
            #include "lexPP.h"
            extern "C" { int atoi(const char*);}
# define YYNEWLINE 10
int lexCode(mn_array *&mn_list){
  int nstr; extern int yyprevious;
  int i, lines=0;
  int ndim=1, ord=0;
  int ps_coef , *base, np=0;
  sht_array *olist=NULL;
  multinomial *mn;
  while((nstr = yylook()) >= 0)
    switch(nstr){
    case 0:
      if(yywrap()) return(0); break;
    case 1:
      {
	lines++;
	if(lines > 1){
	  if(lines%2 == 0){
	    mn =  new multinomial (*olist);
	    ord=0;
	    for(i=0; i<mn->cosize; i++)
	      mn->bp[i]=0.0;}
	  else{
	    ((multinomial **) *mn_list)[np] = mn;
	    np++;}
	}
	
      }
      break;
    case 2:
      {
	if(lines==0){    
	  ndim=atoi(yytext);
	  olist = new sht_array(ndim);
	  base = new int[ndim];
	  mn_list = new mn_array(ndim);}
	
	else if(lines%2 == 1){
	  (*olist)[ord]=atoi(yytext);
	  base[ord]=(*olist)[ord]+1;
	  ord++;}
	else { mn->bp[0] = atoi(yytext);}
	
      }
      break;
    case 3:
      {
	mn->bp[0] = atoi(yytext);
      }
      break;
    case 4:
      {
	
	mn->bp[0] = atof(yytext);
      }
      break;
    case 5:
      {
	mn->bp[0] = atorv();
      }
      break;
    case 6:
      {
	ps_coef = atop(ndim, base);
	if(yytext[0]=='-')
	  mn->bp[ps_coef] = 0.0-1.0;
	else
	  mn->bp[ps_coef] = 1.0;
      }
      break;
    case 7:
      {
	ps_coef = atop(ndim, base);
	mn->bp[ps_coef] = atoi(yytext);
      }
      break;
    case 8:
      {
	ps_coef = atop(ndim, base);
	mn->bp[ps_coef] = atoi(yytext);
      }
      break;
    case 9:
      {
	ps_coef = atop(ndim, base);
	mn->bp[ps_coef] = atof(yytext);
      }
      break;
    case 10:
      {
	ps_coef = atop(ndim, base);
	mn->bp[ps_coef] = atoi(yytext);
      }
      break;
    case 11:
      {
	ps_coef = atop(ndim, base);
	mn->bp[ps_coef] = atof(yytext);
      }
      break;
    case 12:
      {
	ps_coef = atop(ndim, base);
	mn->bp[ps_coef] = atoi(yytext);
      }
      break;
    case 13:
      {
	ps_coef = atop(ndim, base);
	mn->bp[ps_coef] = atof(yytext);
      }
      break;
    case 14:
      {
	ps_coef =  atop(ndim, base);
	if(yytext[0]=='-')
	  mn->bp[ps_coef] = 0.0-1.0;
	else
	  mn->bp[ps_coef] = 1.0;
      }
      break;
    case 15:
      {
	ps_coef = atop(ndim, base);
	if(yytext[0] == 'x')
	  mn->bp[ps_coef] = 1.0;
	else
	  mn->bp[ps_coef] = atoi(yytext);
      }
      break;
    case 16:
      {
	ps_coef = atop(ndim, base);
	if(yytext[0] == 'x')
	  mn->bp[ps_coef] = 1.0;
	else
	  mn->bp[ps_coef] = atoi(yytext);
      }
      break;
case 17:
  {
    ps_coef = atop(ndim, base);
    if(yytext[0] == 'x')
      mn->bp[ps_coef] = 1.0;
    else{
      mn->bp[ps_coef] = atof(yytext);}
  }
  break;
    case 18:
      {
	return 0;
      }
      //break;
    case -1:
      break;
    default:
      fprintf(yyout,"bad switch yylook %d",nstr);
    } 
return(0); 
}
/* end of yylex */

/*  convert the yytext to the value of double */
double  atofv()
{
    int i, j, k, l;
    double a=1.0;
    
    j = atoi(yytext);
    for(i=0; i<yyleng; i++){
	if( yytext[i] == '.' ){
          k = atoi(&yytext[i+1]);
          l = i+1;}}
    for(i=l; i<yyleng+1; i++){
      a *= (double)k*0.1; } 
    if( j < 0 )
      a = (double)j - a;
    else  
      a += (double)j;
    return a;
}

#ifndef USE_INTERVAL
/*  convert the Rat yytext to the value of double */
double  atorv()
{
    int i, j, k;
    double a=0.0;

    j = atoi(yytext);
    for(i=0; i<yyleng; i++){
	if( yytext[i] == '/' ){
          k = atoi(&yytext[i+1]);}}
 
    a = (double)j/(double)k;
    return a;
}
#else
/*  convert the Interval yytext to the value of double */
real  atorv()
{
    int i, j, k;
    real a=0.0;

    j = atoi(yytext);
    for(i=0; i<yyleng; i++){
	if( yytext[i] == '/' ){
          k = atoi(&yytext[i+1]);}}
 
    a = (real)j/(real)k;
    return a;
}
#endif

/*  convert the yytext to proper value of position in bp for polynomial */

int atop(int n, int *base)
{
    int *foot, *power;
    int *ps=new int[n];
    int i, j=0, k, l, ps_bp;

    for(i=0; i<n; i++)
	ps[i]=0;
    for(i=0; i<yyleng; i++){
      if(yytext[i] == 'x' ) j++;}
    foot = new int[j];
    power = new int[j];
    k=0;

    for(i=0; i<yyleng; i++){
      if(yytext[i] == 'x' ){
         foot[k] = atoi(&yytext[i+1]);
         power[k] = 1;                        /* default value of power */
       
	 for(l=i+1; l<yyleng; l++){               /* between the two x's */
           if(yytext[l] == 'x') break;
	   if(yytext[l] == '^'){                  /*  if there is a '^' */
	     power[k] = atoi(&yytext[l+1]);      /* yes, reset power[k] */
	   }}

	 ps[(foot[k])-1]=power[k];
	 k++;
       }
    }
    delete[] foot;
    delete[] power;
    ps_bp = rsort(n, ps, base);
    delete[] ps;
   
    return ps_bp;
}
/* get the position in a string of contpts according to its *ps in base */
/* it is the reverse process of hash                                    */
/* it is: 1+ps[0]+base[0]*ps[1]+base[1]*ps[2]+...+bsae[n-1]*ps[n]       */

int sort(int n, int *ps, int *base)
{

  int i=0; 
  int j=ps[0];

  for(i=1; i<n; i++)
    j +=get_total(i-1, base) * ps[i];

  return j;
}


int rsort(int n, int *ps, int *base)
{
  int i;
  int *rbase=new int[n];
  for(i=0; i<n; i++)
    rbase[n-i-1]=base[i];


  int j=ps[n-1];

  for(i=1; i<n; i++)
    j +=get_total(i, rbase) * ps[n-i-1];

  delete[] rbase;
  return j;
}


/* get the total number of n strings in the base */
/* It is: base[0]*base[1]* .. *base[n-1]         */

int get_total(int n, int *base)
{
  int j=1;

  for(int i=0; i<n; i++)
    j *= base[i];

  return j;
}
int yyvstop[] = {
0,

1,
0,

2,
3,
0,

3,
0,

6,
15,
16,
17,
0,

6,
0,

4,
0,

5,
0,

7,
8,
12,
15,
16,
0,

18,
0,

8,
12,
16,
0,

14,
15,
16,
17,
0,

6,
15,
16,
17,
0,

14,
0,

6,
0,

9,
13,
17,
0,

10,
15,
16,
0,

12,
15,
16,
0,

10,
16,
0,

12,
16,
0,

15,
16,
17,
0,

15,
16,
17,
0,

11,
17,
0,

13,
17,
0,

15,
16,
0,

15,
16,
0,

16,
0,

16,
0,

15,
16,
17,
0,

17,
0,

17,
0,
0};
# define YYTYPE char
struct yywork { YYTYPE verify, advance; } yycrank[] = {
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	1,3,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
1,4,	0,0,	1,4,	0,0,	
0,0,	1,5,	1,5,	1,5,	
1,5,	1,5,	1,5,	1,5,	
1,5,	1,5,	1,5,	4,8,	
4,8,	4,8,	4,8,	4,8,	
4,8,	4,8,	4,8,	4,8,	
4,8,	5,10,	5,11,	5,5,	
5,5,	5,5,	5,5,	5,5,	
5,5,	5,5,	5,5,	5,5,	
5,5,	7,14,	7,14,	7,14,	
7,14,	7,14,	7,14,	7,14,	
7,14,	7,14,	7,14,	8,10,	
8,11,	0,0,	0,0,	0,0,	
14,21,	16,24,	0,0,	19,27,	
23,31,	13,20,	1,6,	30,39,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	6,13,	
0,0,	0,0,	0,0,	0,0,	
35,43,	0,0,	0,0,	0,0,	
0,0,	1,7,	14,22,	16,25,	
17,26,	19,28,	23,32,	29,38,	
33,42,	30,22,	34,25,	4,9,	
9,16,	9,16,	9,16,	9,16,	
9,16,	9,16,	9,16,	9,16,	
9,16,	9,16,	35,44,	5,12,	
10,17,	10,17,	10,17,	10,17,	
10,17,	10,17,	10,17,	10,17,	
10,17,	10,17,	11,18,	11,18,	
11,18,	11,18,	11,18,	11,18,	
11,18,	11,18,	11,18,	11,18,	
36,45,	8,15,	12,19,	12,19,	
12,19,	12,19,	12,19,	12,19,	
12,19,	12,19,	12,19,	12,19,	
15,23,	15,23,	15,23,	15,23,	
15,23,	15,23,	15,23,	15,23,	
15,23,	15,23,	21,29,	21,29,	
21,29,	21,29,	21,29,	21,29,	
21,29,	21,29,	21,29,	21,29,	
22,30,	22,30,	22,30,	22,30,	
22,30,	22,30,	22,30,	22,30,	
22,30,	22,30,	24,33,	24,33,	
24,33,	24,33,	24,33,	24,33,	
24,33,	24,33,	24,33,	24,33,	
25,34,	25,34,	25,34,	25,34,	
25,34,	25,34,	25,34,	25,34,	
25,34,	25,34,	26,35,	26,35,	
26,35,	26,35,	26,35,	26,35,	
26,35,	26,35,	26,35,	26,35,	
27,36,	27,36,	27,36,	27,36,	
27,36,	27,36,	27,36,	27,36,	
27,36,	27,36,	28,37,	28,37,	
28,37,	28,37,	28,37,	28,37,	
28,37,	28,37,	28,37,	28,37,	
31,40,	31,40,	31,40,	31,40,	
31,40,	31,40,	31,40,	31,40,	
31,40,	31,40,	32,41,	32,41,	
32,41,	32,41,	32,41,	32,41,	
32,41,	32,41,	32,41,	32,41,	
37,46,	38,47,	38,47,	38,47,	
38,47,	38,47,	38,47,	38,47,	
38,47,	38,47,	38,47,	39,48,	
39,48,	39,48,	39,48,	39,48,	
39,48,	39,48,	39,48,	39,48,	
39,48,	40,49,	41,50,	48,56,	
51,24,	52,59,	37,28,	42,51,	
42,51,	42,51,	42,51,	42,51,	
42,51,	42,51,	42,51,	42,51,	
42,51,	43,52,	43,52,	43,52,	
43,52,	43,52,	43,52,	43,52,	
43,52,	43,52,	43,52,	55,45,	
41,32,	44,53,	44,53,	44,53,	
44,53,	44,53,	44,53,	44,53,	
44,53,	44,53,	44,53,	45,54,	
45,54,	45,54,	45,54,	45,54,	
45,54,	45,54,	45,54,	45,54,	
45,54,	46,55,	46,55,	46,55,	
46,55,	46,55,	46,55,	46,55,	
46,55,	46,55,	46,55,	47,21,	
49,57,	49,57,	49,57,	49,57,	
49,57,	49,57,	49,57,	49,57,	
49,57,	49,57,	50,58,	50,58,	
50,58,	50,58,	50,58,	50,58,	
50,58,	50,58,	50,58,	50,58,	
53,60,	54,46,	57,50,	58,49,	
63,59,	47,56,	56,61,	56,61,	
56,61,	56,61,	56,61,	56,61,	
56,61,	56,61,	56,61,	56,61,	
59,62,	59,62,	59,62,	59,62,	
59,62,	59,62,	59,62,	59,62,	
59,62,	59,62,	53,44,	54,45,	
57,49,	60,63,	60,63,	60,63,	
60,63,	60,63,	60,63,	60,63,	
60,63,	60,63,	60,63,	61,39,	
62,60,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	61,56,	62,59,	0,0,	
0,0};
struct yysvf yysvec[] = {
0,	0,	0,
yycrank+1,	0,		0,	
yycrank+0,	yysvec+1,	0,	
yycrank+0,	0,		yyvstop+1,
yycrank+11,	0,		0,	
yycrank+23,	0,		yyvstop+3,
yycrank+1,	0,		0,	
yycrank+33,	0,		0,	
yycrank+45,	yysvec+4,	yyvstop+6,
yycrank+84,	0,		0,	
yycrank+96,	0,		0,	
yycrank+106,	0,		0,	
yycrank+118,	0,		0,	
yycrank+1,	0,		0,	
yycrank+2,	yysvec+7,	yyvstop+8,
yycrank+128,	0,		0,	
yycrank+3,	yysvec+9,	yyvstop+13,
yycrank+4,	yysvec+10,	yyvstop+15,
yycrank+0,	yysvec+11,	yyvstop+17,
yycrank+5,	yysvec+12,	yyvstop+19,
yycrank+0,	0,		yyvstop+25,
yycrank+138,	0,		0,	
yycrank+148,	0,		0,	
yycrank+6,	yysvec+15,	yyvstop+27,
yycrank+158,	0,		0,	
yycrank+168,	0,		0,	
yycrank+178,	0,		0,	
yycrank+188,	0,		0,	
yycrank+198,	0,		0,	
yycrank+7,	yysvec+21,	yyvstop+31,
yycrank+9,	yysvec+22,	yyvstop+36,
yycrank+208,	0,		0,	
yycrank+218,	0,		0,	
yycrank+8,	yysvec+24,	yyvstop+41,
yycrank+10,	yysvec+25,	yyvstop+43,
yycrank+22,	yysvec+26,	yyvstop+45,
yycrank+44,	yysvec+27,	yyvstop+49,
yycrank+182,	yysvec+28,	yyvstop+53,
yycrank+229,	0,		0,	
yycrank+239,	0,		0,	
yycrank+177,	yysvec+31,	yyvstop+57,
yycrank+204,	yysvec+32,	yyvstop+60,
yycrank+255,	0,		0,	
yycrank+265,	0,		0,	
yycrank+277,	0,		0,	
yycrank+287,	0,		0,	
yycrank+297,	0,		0,	
yycrank+261,	yysvec+38,	yyvstop+63,
yycrank+179,	yysvec+39,	yyvstop+67,
yycrank+308,	0,		0,	
yycrank+318,	0,		0,	
yycrank+206,	yysvec+42,	0,	
yycrank+181,	yysvec+43,	yyvstop+71,
yycrank+282,	yysvec+44,	yyvstop+74,
yycrank+283,	yysvec+45,	yyvstop+77,
yycrank+203,	yysvec+46,	yyvstop+80,
yycrank+334,	0,		0,	
yycrank+284,	yysvec+49,	yyvstop+83,
yycrank+259,	yysvec+50,	yyvstop+85,
yycrank+344,	0,		0,	
yycrank+357,	0,		0,	
yycrank+321,	yysvec+56,	yyvstop+87,
yycrank+322,	yysvec+59,	yyvstop+91,
yycrank+260,	yysvec+60,	yyvstop+93,
0,	0,	0};
struct yywork *yytop = yycrank+442;
struct yysvf *yybgin = yysvec+1;
char yymatch[] = {
00  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,'+' ,01  ,'+' ,01  ,01  ,
'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,
'0' ,'0' ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
0};
char yyextra[] = {
0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,
0};
/* #ident	"@(#)libl:lib/ncform	1.3" */
#ident	"$Header: /src_trees/lonestar/att/usr/src/lib/libl/RCS/ncform,v 1.5 1991/08/06 11:04:55 jfw Exp $"
int yylineno =1;
# define YYU(x) x
# define NLSTATE yyprevious=YYNEWLINE
char yytext[YYLMAX];
struct yysvf *yylstate [YYLMAX], **yylsp, **yyolsp;
char yysbuf[YYLMAX];
char *yysptr = yysbuf;
int *yyfnd;
extern struct yysvf *yyestate;
int yyprevious = YYNEWLINE;
int yyback(int * p, int m) {
if (p==0) return(0);
while (*p)
	{
	if (*p++ == m)
		return(1);
	}
return(0);
}
	/* the following are only used in the lex library */
int yyinput(){
	return(input());
	}
void yyoutput(int c) {
	output(c);
	}
void yyunput(int c) {
	unput(c);
	}

int yylook(){
	register struct yysvf *yystate, **lsp;
	register struct yywork *yyt;
	struct yysvf *yyz;
	int yych, yyfirst;
	struct yywork *yyr;
# ifdef LEXDEBUG
	int debug;
# endif
	char *yylastch;
	/* start off machines */
# ifdef LEXDEBUG
	debug = 0;
# endif
	yyfirst=1;
	if (!yymorfg)
		yylastch = yytext;
	else {
		yymorfg=0;
		yylastch = yytext+yyleng;
		}
	for(;;){
		lsp = yylstate;
		yyestate = yystate = yybgin;
		if (yyprevious==YYNEWLINE) yystate++;
		for (;;){
# ifdef LEXDEBUG
			if(debug)fprintf(yyout,"state %d\n",yystate-yysvec-1);
# endif
			yyt = yystate->yystoff;
			if(yyt == yycrank && !yyfirst){  /* may not be any transitions */
				yyz = yystate->yyother;
				if(yyz == 0)break;
				if(yyz->yystoff == yycrank)break;
				}
			*yylastch++ = yych = input();
			yyfirst=0;
		tryagain:
# ifdef LEXDEBUG
			if(debug){
				fprintf(yyout,"char ");
				allprint(yych);
				putchar('\n');
				}
# endif
			yyr = yyt;
			if ( (int)yyt > (int)yycrank){
				yyt = yyr + yych;
				if (yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transitions */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					goto contin;
					}
				}
# ifdef YYOPTIM
			else if((int)yyt < (int)yycrank) {		/* r < yycrank */
				yyt = yyr = yycrank+(yycrank-yyt);
# ifdef LEXDEBUG
				if(debug)fprintf(yyout,"compressed state\n");
# endif
				yyt = yyt + yych;
				if(yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transitions */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					goto contin;
					}
				yyt = yyr + YYU(yymatch[yych]);
# ifdef LEXDEBUG
				if(debug){
					fprintf(yyout,"try fall back character ");
					allprint(YYU(yymatch[yych]));
					putchar('\n');
					}
# endif
				if(yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transition */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					goto contin;
					}
				}
			if ((yystate = yystate->yyother) && (yyt= yystate->yystoff) != yycrank){
# ifdef LEXDEBUG
				if(debug)fprintf(yyout,"fall back to state %d\n",yystate-yysvec-1);
# endif
				goto tryagain;
				}
# endif
			else
				{unput(*--yylastch);break;}
		contin:
# ifdef LEXDEBUG
			if(debug){
				fprintf(yyout,"state %d char ",yystate-yysvec-1);
				allprint(yych);
				putchar('\n');
				}
# endif
			;
			}
# ifdef LEXDEBUG
		if(debug){
			fprintf(yyout,"stopped at %d with ",*(lsp-1)-yysvec-1);
			allprint(yych);
			putchar('\n');
			}
# endif
		while (lsp-- > yylstate){
			*yylastch-- = 0;
			if (*lsp != 0 && (yyfnd= (*lsp)->yystops) && *yyfnd > 0){
				yyolsp = lsp;
				if(yyextra[*yyfnd]){		/* must backup */
					while(yyback((*lsp)->yystops,-*yyfnd) != 1 && lsp > yylstate){
						lsp--;
						unput(*yylastch--);
						}
					}
				yyprevious = YYU(*yylastch);
				yylsp = lsp;
				yyleng = yylastch-yytext+1;
				yytext[yyleng] = 0;
# ifdef LEXDEBUG
				if(debug){
					fprintf(yyout,"\nmatch ");
					sprint(yytext);
					fprintf(yyout," action %d\n",*yyfnd);
					}
# endif
				return(*yyfnd++);
				}
			unput(*yylastch);
			}
		if (yytext[0] == 0  /* && feof(yyin) */)
			{
			yysptr=yysbuf;
			return(0);
			}
		yyprevious = yytext[0] = input();
		if (yyprevious>0)
			output(yyprevious);
		yylastch=yytext;
# ifdef LEXDEBUG
		if(debug)putchar('\n');
# endif
		}
	}


