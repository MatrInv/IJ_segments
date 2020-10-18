import java.awt.*;

public class DigitalLineSegment {
	// droite 4-connexe du 1er quadrant (déplacements 0 et 1)
	// mu <= ax - by < mu + a+ b 

	private int _a, _b, _mu; /* caracteristiques */
	private Point _U; /*pt d'appui sup d'absc min*/
	private Point _L; /*pt d'appui inf d'absc min*/
	private Point _Us;/*pt d'appui sup d'absc max*/
	private Point _Ls;/*pt d'appui inf d'absc max*/


	DigitalLineSegment(int firstMove){
		_U = new Point();
		_L = new Point();

		_mu = 0;

		if (firstMove == 0){
			_a = 0;
			_b = 1;
			_Us = new Point(1, 0);
			_Ls = new Point(1, 0);
		}
		else{
			_a = 1;
			_b = 0;
			_Us = new Point(0, 1);
			_Ls = new Point(0, 1);
		}
	}

	DigitalLineSegment(DigitalLineSegment d){
		_a = d._a;
		_b = d._b;
		_mu = d._mu;
		_U = new Point(d._U);
		_L = new Point(d._L);
		_Us = new Point(d._Us);
		_Ls = new Point(d._Ls);
	}

	DigitalLineSegment(){
		_a = 0;
		_b = 0;
		_mu = 0;
		_U = new Point();
		_L = new Point();
		_Us = new Point();
		_Ls = new Point();
	}


	public void copy(DigitalLineSegment d){
		_a = d._a;
		_b = d._b;
		_mu = d._mu;
		_U.setLocation(d._U);
		_L.setLocation(d._L);
		_Us.setLocation(d._Us);
		_Ls.setLocation(d._Ls);
	}

	public int getA(){
		return _a;
	}
	public int getB(){
		return _b;
	}

	

	// Ajout d'un point a l'extremite d'abscisse positive d'un segment
	// de droite 4-connexe du 1er quadrant 

	boolean addPositivePoint (int x, int y){
		int r;
		r = _a*x-_b*y;
		if(_mu<=r && r<_mu+_a+_b){
			if (r==_mu)
				_Us.move(x, y);
			if(r==_mu+_a+_b-1)
				_Ls.move(x, y);
			return true;
		}

		if(r==_mu-1){ // la pente augmente
			//mise a jour des points d'appui
			_L.setLocation(_Ls);
			_Us.move(x, y);

			//mise a jour des caracteristiques
			_a = y - (int)_U.getY();
			_b = x - (int)_U.getX();
			_mu = _a*x-_b*y;
			return true;
		}

		if(r==_mu+_a+_b){ // la pente diminue
			//mise a jour des points d'appui
			_U.setLocation(_Us);
			_Ls.move(x, y);

			//mise a jour des caracteristiques
			_a = y - (int)_L.getY();
			_b = x - (int)_L.getX();
			_mu = _a*(x-1)-_b*(y+1)+1;
			return true;
		}
		// sinon le nouveau point ne permet pas de continuer la droite
		return false;
	}


	// Ajout d'un point a l'extremite d'abscisse negative d'un segment
	// de droite 4-connexe du 1er quadrant 

	boolean addNegativePoint(int x, int y){
		int r;
		r = _a*x-_b*y;
		if(_mu<=r && r< _mu+_a+_b){
			if (r==_mu)
				_U.move(x, y);
			if(r == _mu+_a+_b-1)
				_L.move(x, y);
			return true;
		}

		if(r==_mu-1){ //la pente diminue
			// mise a jour des points d'appui
			_Ls.setLocation(_L);
			_U.move(x, y);

			// mise a jour des caracteristiques
			_a = (int)_Us.getY() - y;
			_b = (int)_Us.getX() - x;
			_mu = _a*x-_b*y;
			return true;
		}

		if(r==_mu+_a+_b){ // la pente augmente
			// mise a jour des points d'appui
			_Us.setLocation(_U);
			_L.move(x, y);

			// mise a jour des caracteristiques
			_a = (int)_Ls.getY()-y;
			_b = (int)_Ls.getX()-x;
			_mu = _a*(x-1)-_b*(y+1)+1;
			return true;
		}
		// le nouveau point ne permet pas de continuer la droite
		return false;
	}

}
