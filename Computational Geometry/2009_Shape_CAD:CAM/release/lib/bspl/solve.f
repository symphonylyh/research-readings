c
c  Copyright (C) 1990 Massachusetts Institute of Technology all rights reserved
c
c
c	........FUNCTION SIMUL FOR SOLVING A SET OF SIMULTANEOUS EQUATIONS .....
c     iwks and y is workspace size
	function simul(n,a,x,eps,indic,nrc,iwks,y)
	implicit real*8(A-H, O-Z)
	real*8 a,x,eps,y,simul
	dimension iwks(3,nrc),y(nrc)
	dimension x(n),a(nrc,nrc)

c*****************************************************************
c
	max = n
	if(indic.ge.0) max = n+1
c
c	......IS N LARGER THAN 50 No limit anymore.......
c
c	if(n.le.50) go to 5
c	write(*,200) 
c	simul = 0.0
c	return
c
c	......BEGIN ELIMINATION PROCEDURE.........
c
5	deter = 1.0
	do 18 k=1,n
	km1 = k-1
c
c	......SEARCH FOR THE PIVOT ELEMENT........
c
	pivot = 0.0
	do 11 i=1,n
	do 11 j=1,n
c
c 	......SCAN IROW AND JCOL ARRAYS FOR INVALID PIVOT SUBSCRIPTS.....
c	
	if(k.eq.1) go to 9
	do 8 iscan = 1, km1
	do 8 jscan = 1, km1
	if(i.eq.iwks(1,iscan)) go to 11
	if(j.eq.iwks(2,jscan)) go to 11
8 	continue
9 	if(dabs(a(i,j)).le.dabs(pivot)) go to 11
	pivot = a(i,j)
	iwks(1,k) = i
	iwks(2,k) = j
11 	continue
c	
c	......INSURE THAT SELECTED PIVOT IS LARGER THAN EPS .........
c	
	if( dabs(pivot) .gt. eps ) go to 13
	simul = 0.0
	return
c
c	......UPDATE THE DETERMINANT VALUE ...........
c
13	irowk = iwks(1,k)
	jcolk = iwks(2,k)
	deter = deter*pivot
c
c	......NORMALIZE PIVOT ROW ELEMENTS ...........
c
	do 14 j=1,max
14 	a(irowk,j) = a(irowk,j)/pivot
c
c	......CARRY OUT ELIMINATION AND DEVELOP INVERSE .........
c
	a(irowk,jcolk) = 1.0/pivot
	do 18 i = 1,n
	aijck = a(i,jcolk)
	if( i.eq.irowk) go to 18
	a(i,jcolk) = - aijck/pivot
	do 17 j=1,max
17	if(j.ne.jcolk) a(i,j) = a(i,j) - aijck*a(irowk,j)
18	continue
c
c	....... ORDER SOLUTION VALUES (IF ANY) AND CREATE JORD ARRAY .....
c
	do 20 i=1,n
	irowi = iwks(1,i)
	jcoli = iwks(2,i)
	iwks(3,irowi) = jcoli
20	if( indic .ge. 0 ) x(jcoli) = a(irowi,max)
c
c	....... ADJUST SIGN OF DETERMINANT .........
c
	intch = 0
	nm1 = n-1
	do 22 i=1,nm1
	ip1 = i+1
	do 22 j=ip1,n
	if( iwks(3,j) .ge. iwks(3,i)) go to 22
	jtemp = iwks(3,j)
	iwks(3,j) = iwks(3,i)
	iwks(3,i) = jtemp
	intch = intch + 1
22 	continue
	if( intch/2*2 .ne. intch ) deter = - deter
c
c	....... IF INDIC IS POSITIVE RETURN WITH RESULTS ...........
c
	if (indic .le. 0) go to 26
	simul = deter
	return
c
c	....... IF INDIC IS NEGATIVE OR ZERO, UNSCRAMBLE THE INVERSE 
c		FIRST BY ROWS .........
c
26	do 28 j=1,n
	do 27 i=1,n
	irowi = iwks(1,i)
	jcoli = iwks(2,i)
27	y(jcoli) = a(irowi,j)
	do 28 i=1,n
28	a(i,j) = y(i)
c
c	....... THEN BY COLUMNS ....... 
c
	do 30 i=1,n
	do 29 j=1,n
	irowj = iwks(1,j)
	jcolj = iwks(2,j)
29 	y(irowj) = a(i,jcolj)
	do 30 j=1,n
30 	a(i,j) = y(j)
c
c	....... RETURN FOR INDIC NEGATIVE OR ZERO ..........
c
	simul = deter
	return
c
c	....... FORMAT FOR OUTPUT STATEMENT  ...........
c
200	format(10h0N too big )
c
	end


