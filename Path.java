import java.awt.*;
/**
*  Classe permettant des calculs d'informations géométriques 
*  en chaque  point discret d’un chemin 4-connexe issu d'un code de Freeman
*/ 
public class Path {

	private class ContourPoint{
		// repère x "vers le haut" et y "vers la droite"
		// (x,y) désigne le point demi-entier (x+1/2, y+1/2)
		int _x;
		int _y; 
		
		// Informations géométriques calculées :
		
		// vecteur tangent symétrique
		int _tx;
		int _ty;
		
		// vecteur tangent  dans la direction positive
		int _tpx;
		int _tpy;

		int _lgrTp; // Nombre de points de la tangente positive
		
		// courbure : TODO
	}

	private int _n;       // nb de points de contour
	private ContourPoint[] _tabPts; // tableau des points de contour
	private int _xm; // abscisse minimum
	private int _ym; // ordonnée minimum
	private int _xM; // abscisse maximum
	private int _yM; // ordonnée maximum


	Path(FreemanCode c){
		_n = c.getLength();
		_tabPts = new ContourPoint[_n];
		for(int i = 0; i<_n; i++)
			_tabPts[i] = new ContourPoint();
		computeCoordinates(c);  // calcul des coordonnées des points du contour c 
		computeTangents(c);		// calcul des tangentes symétriques en chaque point
		computeTangentsPos(c);  // calcul des tangentes positives en chaque point
	}

	/*
	* Calcul les points du contour à partir du code de Freeman c
	*/
	private void computeCoordinates(FreemanCode c){
		int x=(int)c.getP0().getX(), y=(int)c.getP0().getY();
		_xm=_xM=x;
		_ym=_yM=y;
		_tabPts[0]._x=x;
		_tabPts[0]._y=y;

		for (int i = 0; i < _n-1; i++){
			x += c.getDeltaX(i);
			if (x < _xm) 
				_xm = x;
			if (x > _xM)
				_xM = x;
			y += c.getDeltaY(i);
			if (y < _ym)
				_ym = y;
			if (y > _yM)
				_yM = y;
			_tabPts[i+1]._x=x;
			_tabPts[i+1]._y=y;
		}
	}
	
	/*
	* Calcul de la tangente symétrique en chaque point du contour
	*/
	private void computeTangents(FreemanCode c){
		DigitalTangent t;
		int[][] sym4 = {{1,1},{-1,1},{-1,-1},{1,-1}};

		for(int i = 0; i < _n; i++){
			t= new DigitalTangent(c, i);
			// calcul de symetries a partir des caracteristiques 
			// rapportees au 1er quadrant
			int a = t.getA();
			int b = t.getB();
			int quad = t.getQuadrant();
			_tabPts[i]._tx = sym4[quad][0]*b;
			_tabPts[i]._ty = sym4[quad][1]*a;
		}

	}

	/*
	* Calcul de la tangente positive en chaque point du contour
	*/
	private void computeTangentsPos(FreemanCode c){
		DigitalTangentPositive tp;
		int[][] sym4 = {{1,1},{-1,1},{-1,-1},{1,-1}};

		for(int i = 0; i < _n; i++){
			tp= new DigitalTangentPositive(c, i);
			// calcul de symetries a partir des caracteristiques 
			// rapportees au 1er quadrant
			int a = tp.getA();
			int b = tp.getB();
			int quad = tp.getQuadrant();
			_tabPts[i]._tpx = sym4[quad][0]*b;
			_tabPts[i]._tpy = sym4[quad][1]*a;
			_tabPts[i]._lgrTp=tp.getN();
		}

	}

	public int getLength(){
		return _n;
	}


	public Point getPMin(){
		return new Point(_xm, _ym);
	}
	public Point getPMax(){
		return new Point(_xM, _yM);
	}

	public int getContourX(int i){
		return _tabPts[i]._x;
	}

	public int getContourY(int i){
		return _tabPts[i]._y;
	}
	public int getContourTX(int i){
		return _tabPts[i]._tx;
	}

	public int getContourTY(int i){
		return _tabPts[i]._ty;
	}
	
	public int getContourTPX(int i){
		return _tabPts[i]._tpx;
	}

	public int getContourTPY(int i){
		return _tabPts[i]._tpy;
	}
	
	public int getContourTPLgr(int i){
		return _tabPts[i]._lgrTp;
	}

}

