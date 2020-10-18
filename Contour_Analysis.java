import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;
import java.util.ArrayList;

public class Contour_Analysis implements PlugInFilter {
    ImagePlus imp;
    int W, H; // taille de l'image initiale
    int iMin, iMax, jMin, jMax; // coordonnées de la ROI contenant l'objet analysé
    int GL; // niveau de gris de l'objet analysé
    int ZL; // niveau de zoom de la visualisation
    private static final int NSIZE = 5;
    private static final double ACCURACY = 100; //variable donnant la précision du tracé des normales


    /**
     * Set ups this filter. Process only 8 bits grey-level images.
     */
    public int setup(String arg, ImagePlus imp) {
	if ( arg == "about" ) {
	    IJ.showMessage( "about Contour_Analysis ...", 
			    "displays the contour's dominant points");
	    return DONE;
	}
	this.imp = imp;
	return DOES_8G;
    }



    private int getIVisu(Path p, int ind){
	return (p.getContourX(ind)-iMin+1)*ZL;
    }

    private int getJVisu(Path p, int ind){
	return (H-1-p.getContourY(ind)-jMin)*ZL;
    }

    private int getIVisu(int x){
	return (x-iMin+1)*ZL;
    }

    private int getJVisu(int y){
	return (H-1-y-jMin)*ZL;
    }

    //fonction triangle/lambda
    private static double triangleFunc(double x){
	double pic = 0.5;

	if(x>=1 || x<=0 )
	    return 0;
	else
	    return (1.-Math.abs(2.*(x-pic)));
    }

    public void run( ImageProcessor ip ) {
	// Permet de recuperer des parametres.
	GenericDialog gd = new GenericDialog( "Contour Extraction parameters",
					      IJ.getInstance() );
	gd.addNumericField( "gray level", 0, 0, 3, "" );
	String[] ZLChoice = {"2", "4", "8", "16"};
	gd.addChoice("zoom level", ZLChoice, "4"); 
	gd.showDialog();
	if ( gd.wasCanceled() ) {
	    IJ.error( "PlugIn cancelled" );
	    return;
	}
	GL = (int) gd.getNextNumber();
	ZL = (int)Math.pow(2, gd.getNextChoiceIndex()+1);
	W = ip.getWidth();
	H = ip.getHeight();

	// Crée le contour de la région
	int margin = 3;
	FreemanCode c = new FreemanCode(ip, GL);
	Path p = new Path(c);
	iMin = (int)p.getPMin().getX() - margin + 1;
	iMax = (int)p.getPMax().getX() + margin;
	jMax = H-1 - ((int)p.getPMin().getY() - margin + 1);
	jMin = H-1 - ((int)p.getPMax().getY() + margin);


	// Crée l'image de sortie
	ip.setRoi(iMin, jMin, iMax-iMin+1, jMax-jMin+1);
	ImageProcessor ipz = ip.resize((iMax-iMin+1)*ZL, (jMax-jMin+1)*ZL);
	ImageProcessor ipz2 = ip.resize((iMax-iMin+1)*ZL, (jMax-jMin+1)*ZL);
	
	// Tracé de la grille 1/2
	ipz.setColor(Color.black);
	for(int i=0; i<W*ZL; i+=ZL) 
	    for(int j=0; j<H*ZL; j+=ZL)
		ipz.drawPixel(i, j);

	// Tracé du chemin
	ImageProcessor ipzc = ipz.convertToRGB();
	ImageProcessor ipzc2 = ipz2.convertToRGB();
	ipzc.setColor(Color.red);
	int ind, i1, j1, i2=0, j2=0;
	for(ind=0; ind<p.getLength()-1; ind++){
	    i1 = getIVisu(p, ind); 
	    j1 = getJVisu(p, ind); 
	    i2 = getIVisu(p, ind+1); 
	    j2 = getJVisu(p, ind+1); 
	    ipzc.drawLine(i1, j1, i2, j2);
	}
	i1 = i2;
	j1 = j2;
	i2 = getIVisu(p, 0); 
	j2 = getJVisu(p, 0); 
	ipzc.drawLine(i1, j1, i2, j2);

	// Tracé des normales associées aux tangentes symétriques
	/*ipzc.setColor(Color.green);
	for(ind = 0; ind < p.getLength(); ind+=5){
	    int x = p.getContourX(ind);
	    int y = p.getContourY(ind);
	    int a = p.getContourTY(ind);
	    int b = p.getContourTX(ind);
	    double norme = Math.sqrt(a*a+b*b);
	    i1 = getIVisu(x);
	    j1 = getJVisu(y);
	    i2 = getIVisu((int)(x-(NSIZE*a)/norme+0.5));
	    j2 = getJVisu((int)(y+(NSIZE*b)/norme+0.5));
	    ipzc.drawLine(i1, j1, i2, j2);

	    }*/

	/************************************************************************
	*CALCUL DE LA COUVERTURE TANGENTIELLE DU CONTOUR
	**************************************************************************/

	int n = c.getLength();
	int seg = 0;
	
	//on garde le premier segment en mémoire
	int end = p.getContourTPLgr(0);

	//liste des indices des points de départ des segments maximaux
	//elle représente la couverture tangentielle
	ArrayList<Integer> tanCoverStart = new ArrayList<Integer>();

	//si le premier segment n'est pas inclu dans le dernier, on ajoute 0 à la couverture tangentielle
	if(end%n != (n-1+p.getContourTPLgr(n-1))%n){
	    tanCoverStart.add(0);
	    IJ.log("MSS_"+seg+" : "+0+" --> "+p.getContourTPLgr(0));
	    seg++;
	}
	
	for(int i = 0; i < n; i++){
	    
	    //si le point i appartient à un autre segment
	    if(i+p.getContourTPLgr(i) != end){
		
		//on garde en mémoire la fin du nouveau segment			     
		end=i+p.getContourTPLgr(i);

		//on range l'indice du point de départ du nouveau segment dans une liste qui représente notre couverture tangentielle
		tanCoverStart.add(i);
		
		IJ.log("MSS_"+ seg +" : "+i+" --> "+end%n+" : "+p.getContourTX(i)+", "+p.getContourTY(i));
		
		seg++;//incrémente le num de segment
	    }			
		
	}
	
	
	/************************************************************************
	*CALCUL DE LA DIRECTION TANGENTIELLE EN CHAQUE POINT 
	*AVEC L'ESTIMATEUR LAMBDA-MST 
	**************************************************************************/
	double dirTanNumerator;
	double dirTanDenominator;
	double[] dirTan = new double[n];
	double excentricity;
	double angle;
	int s;
	
	//Pour tout i<n
	for(int i=0;i<n;i++){
	    dirTanNumerator = 0;
	    dirTanDenominator = 0;
	    //et pour tout segment s
	    for(int si=0;si<tanCoverStart.size();si++){

		s = tanCoverStart.get(si);

		//si le point Ai appartient au segment s
		if( (i>=s && i<= s+p.getContourTPLgr(s)) || (s+p.getContourTPLgr(s)>=n && (s+p.getContourTPLgr(s))%n>=i) ){
		    
		    //on calcule l'excentricité de s au point Ai
		    //excentricity = triangleFunc( (double)(Math.abs(p.getContourX(i)-p.getContourX(s)) + Math.abs(p.getContourY(i)-p.getContourY(s))) / (double)(Math.abs(p.getContourX((s+p.getContourTPLgr(s))%n)-p.getContourX(s)) + Math.abs(p.getContourY((s+p.getContourTPLgr(s))%n)-p.getContourY(s))) );
		    excentricity = triangleFunc( (double)(i-s) / (double)(p.getContourTPLgr(s)) );
		    
		    //et l'angle du segment
		    angle = Math.atan2((double)(p.getContourY((s+p.getContourTPLgr(s))%n)-p.getContourY(s)), (double)(p.getContourX((s+p.getContourTPLgr(s))%n)-p.getContourX(s))) ;

		    //on met à jour les sommes permettant de calculer la direction tangentielle
		    dirTanNumerator += (excentricity * angle);
		    dirTanDenominator += excentricity;
			
		}
			
	    }

	    if(dirTanDenominator != 0){
		dirTan[i] = dirTanNumerator/dirTanDenominator;
		IJ.log("dirTan "+i+" : "+dirTan[i]);
	    }
	    else
		IJ.log("Division par 0 impossible : vous n'avez pas de contour !");
	}
	

	/************************************************************************
	*POLYGONALISATION DU CONTOUR
	**************************************************************************/

	ArrayList<Integer> polygon = new ArrayList<Integer>(); //on représente le polygone comme une liste d'indice de point dans le contour
	int current = 0;
	while(current<n){
	    current+=p.getContourTPLgr(current);
	    if(current<n)
		polygon.add(current);
	    else if(current%n != p.getContourTPLgr(0))
		polygon.add(0);     
	}

	double perimetre = 0; //le perimetre du polygone
	int diffX, diffY;
	
	for(int i=0;i<polygon.size();i++){
	    diffX = (p.getContourX(i)-p.getContourX((i+1)%polygon.size()));
	    diffY = (p.getContourY(i)-p.getContourY((i+1)%polygon.size()));

	    perimetre += Math.sqrt(diffX*diffX + diffY*diffY);
	}
	
	
	//*********************************************************************************	
	// Tracé des normales associées à l'estimateur lambda-mst
	double diffAngle;
	for(ind = 0; ind < p.getLength(); ind+=1){
	    
	    diffAngle = (((dirTan[ind]+2*Math.PI)%(2*Math.PI))-((Math.atan2((double)p.getContourTPY(ind),(double)p.getContourTPX(ind))+2*Math.PI)%(2*Math.PI)));

	    /*if(ind == 98 || ind == 151|| ind==333 || ind ==633 || ind ==800){
		ipzc.setColor(Color.red);
		IJ.log(ind+" : LambdaMST="+(dirTan[ind]+Math.PI/2.)+" TPositive="+ (Math.atan2((double)p.getContourTPY(ind),(double)p.getContourTPX(ind))+Math.PI/2.) );
	    }
	    else*/ if(Math.abs(diffAngle)>Math.PI*1./6. && Math.abs(diffAngle)<Math.PI*11./6.){
		ipzc.setColor(Color.blue);
	    }
	    else
	    ipzc.setColor(Color.green);
	    
	    int x = p.getContourX(ind);
	    int y = p.getContourY(ind);
	    int a = (int)Math.round( Math.cos( dirTan[ind] +Math.PI/2. ) * ACCURACY );
	    int b = (int)Math.round( Math.sin( dirTan[ind] +Math.PI/2. ) * ACCURACY );
	    double norme = Math.sqrt(a*a+b*b);
	    i1 = getIVisu(x);
	    j1 = getJVisu(y);
	    i2 = getIVisu((int)(x+(NSIZE*a)/norme+0.5));
	    j2 = getJVisu((int)(y+(NSIZE*b)/norme+0.5));
	    ipzc.drawLine(i1, j1, i2, j2);

	}

	//*********************************************************************************	
	// Tracé du polygone
	ipzc2.setColor(Color.red);

	for(ind = 0; ind < polygon.size(); ind++){
	    
	    int x = p.getContourX(polygon.get(ind));
	    int y = p.getContourY(polygon.get(ind));
	    int a = p.getContourX(polygon.get((ind+1)% polygon.size()));
	    int b = p.getContourY(polygon.get((ind+1)% polygon.size()));
	    i1 = getIVisu(x);
	    j1 = getJVisu(y);
	    i2 = getIVisu(a);
	    j2 = getJVisu(b);
	    ipzc2.drawLine(i1, j1, i2, j2);

	}
	IJ.log("Perimetre = "+perimetre);
	
	//*********************************************************************************
	
	ImagePlus visu = new ImagePlus("Discrete contour + normals", ipzc);
	ImagePlus visu2 = new ImagePlus("Polygonized", ipzc2);
	
	// displays it in a window.
	visu.show();
	visu2.show();
	// forces the redisplay.
	visu.updateAndDraw();
	visu2.updateAndDraw();
    }

}
			      
