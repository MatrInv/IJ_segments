public class DigitalTangent {

	private DigitalLineSegment _d; // droite 4-c du premier quadrant 
	private int _n; // tangente symetrique de taille 2n+1
	private int _dir0, _dir1; // directions paire et impaire composant la tangente, 
	// l'une d'elle peut etre indefinie = -1
	private int _xm, _ym, _xM, _yM; // points extremes de la tangente dans le 1er quadrant
	// axes x vers la droite et y vers le haut


	private boolean consistentMove (int dir) { // remplit si nécessaire _dir0 ou _dir1
		if (dir == _dir0 || dir == _dir1)
			return true;

		if (_dir0 == -1 && (dir == _dir1-1 || dir == (_dir1+1)%4)){
			_dir0 = dir;
			return true;
		}

		if (_dir1 == -1 && (dir == _dir0+1 || dir == (_dir0+3)%4)){
			_dir1 = dir;
			return true;
		}
		return false;
	}


	// calcul de tangente symétrique par ajout de paires de points
	DigitalTangent(FreemanCode c, int ind){
		int dir;
		int xM, yM, xm, ym;

		// initialisations
		_n = 0;
		dir = c.getCode(ind, 1);

		if (dir%2 == 0){ // pente initiale nulle
			xM = 1; 
			yM = 0;
			_dir0 = dir;
			_dir1 = -1;
			_d = new DigitalLineSegment(0);
		}
		else{ // pente initiale infinie
			xM = 0;
			yM = 1;
			_dir0 = -1;
			_dir1 = dir;
			_d = new DigitalLineSegment(1);
		}

		dir = c.getCode(ind, -1);
		consistentMove(dir); // tjrs vrai
		DigitalLineSegment dSave = new DigitalLineSegment();

		if(dir%2 == 0){
			xm = -1; 
			ym = 0;
		}
		else{
			xm = 0; 
			ym = -1;
		}

		boolean ok = _d.addNegativePoint(xm, ym);
		int i = 2;
		while(ok){
			dSave.copy(_d);
			_xm = xm; _ym = ym;
			_xM = xM; _yM = yM;
			_n++;

			// Ajout du point positif 
			dir = c.getCode(ind, i);
			if (consistentMove(dir)){
				xM += dir%2 == 0 ? 1 : 0;
				yM += dir%2 == 1 ? 1 : 0;
				ok = _d.addPositivePoint(xM, yM);
			}
			else
				ok = false;

			// Ajout du point negatif */
			dir = c.getCode(ind, -i);
			if (ok && consistentMove(dir)){
				xm -= dir%2 ==0 ? 1 : 0;
				ym -= dir%2 ==1 ? 1 : 0;
				ok = _d.addNegativePoint(xm, ym);
			}
			else 
				ok = false;

			i++;
		}
		_d.copy(dSave);
	}


	public int getA(){
		return _d.getA();
	}

	public int getB(){
		return _d.getB();
	}

	public int getQuadrant(){
		if (_dir0 == -1)
			return _dir1;
		if (_dir1 == -1)
			return _dir0;
		if (_dir0 == 0 && _dir1 == 3)
			return 3;
		return _dir0 <= _dir1 ? _dir0 : _dir1;
	}

}
