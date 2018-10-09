C     file di variabili common utili per il bario e stronzio e lantanio e zinco!  

      	integer count
	real zbario(200),ba1(200),ba2(200)
        common/bar/zbario,ba1,ba2
	common/bar2/count

  	integer countla
	real zlanta(200),la1(200),la2(200)
        common/lan/zlanta,la1,la2
	common/lan/countla
  

	real mass(100),zetael(100),yield(6,6)
        integer row,col
        common/sr1/row,col
        common/sr2/mass,zetael,yield

	real masszn(100),zetaelzn(100),yieldzn(60,60)
        integer rowzn,colzn
        common/zn1/rowzn,colzn
        common/zn2/masszn,zetaelzn,yieldzn
	
	real massfe(100),zetaelfe(100),yieldfe(60,60)
        integer rowfe,colfe
        common/fe1/rowfe,colfe
        common/fe2/massfe,zetaelfe,yieldfe

	real yieldmn(60,60)
        common/mn/yieldmn