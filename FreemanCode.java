import java.util.*; 
import ij.process.*;
import java.awt.*;

public class FreemanCode {
    
    private ArrayList<Integer> _code;
    private int _n;  // nb de directions
    private Point _p0; // premier point obtenu par balayage de l'image haut-bas/gauche-droite

    private static final int deltaI[] = {1, 0, -1, 0};
    private static final int deltaJ[] = {0, -1, 0, 1};
    private static final int prev[] = {3, 0, 1, 2};
    private static final int next[] = {1, 2, 3, 0};
    private static final int initialCapacity = 100;

    private static final int deltaX[] = {1, 0, -1, 0};
    private static final int deltaY[] = {0, 1, 0, -1};


    //Construction d'un code de Freeman par suivi de la composante connexe de niveau 
    //de gris ng dans l'image (composante 8-connexe). Hypothèse : bordure entourant l'image
    FreemanCode(ImageProcessor ip, int gl){
	_code = new ArrayList<Integer>(initialCapacity);
	_p0 = new Point();
	_n = 0;
	int i, j, i0=0, j0=0;
	int k = 0;
	int w = ip.getWidth();
	int h = ip.getHeight();
	boolean found = false;
	// parcours de l'image jusqu'à trouver le premier pixel de la composante connexe
	for(j=0; j<h && !found; j++)
	    for(i=0; i<w && !found; i++){
		if (ip.getPixel(i, j) == gl){
		    found = true;
		    i0 = i;
		    j0 = j;
		}
	    }
	// premières directions
	_code.add(k++, 2);
	_code.add(k++, 3);

	i = i0;
	j = j0;
	int dir = 3;

	while (i != i0 || j != j0 || dir != 2){ // on n'est pas revenu au 1er point
	    int i_test = i+deltaI[dir]+deltaI[prev[dir]];
	    int j_test = j+deltaJ[dir]+deltaJ[prev[dir]];
	    if(ip.getPixel(i_test, j_test) == gl){
		dir = prev[dir];
		i = i_test;
		j = j_test;
	    }	
	    
	    else{
		i_test = i+deltaI[dir];
		j_test = j+deltaJ[dir];
		if(ip.getPixel(i_test, j_test) == gl){
		    i = i_test;
		    j = j_test;
		}	
		else
		    dir	= next[dir];
	    }
	    _code.add(k++,dir);
	}
	_code.remove(k-1);
	_n = k - 1; //on est revenu au 1er pas
	_p0.move(i0, ip.getHeight()-1-j0);
    }	

    // taille du chemin en nombre de directions élémentaires
    int getLength (){
	return _n;
    }
    // premier point
    Point getP0 (){
	return new Point(_p0);
    }

    // Variations en x et y pour un pas du chemin defini par le code 
    // l'indice ind doit vérifier 0 <= i < _n
    // Repère utilisé
    //
    //   ^ y
    //   |
    //   |
    //    ------> x 

    int getDeltaX (int ind){
	return deltaX[_code.get(ind)];
    }
    int getDeltaY (int ind){
	return deltaY[_code.get(ind)];
    }

    // direction élémentaire - hypothèse : chemin fermé
    int getCode (int ind, int shift){
	if(shift > 0)
	    return _code.get((ind+shift-1)%_n);
	if(ind+shift >= 0)
	    return _code.get(ind+shift);
	return _code.get((ind+shift)%_n+_n);
    }
}
