package algorithms;

import java.awt.Point;

public class Tools {
	
	// Droite perpendiculaire a line et passant par son milieu
	public MyDroite droitePerpendiculairePassantParMilieu(Point p, Point q) {
		double a;
		Point milieu;
		milieu = new Point((int) ((p.getX() + q.getX()) / 2), (int)((p.getY() + q.getY()) / 2));
		a = (p.getY() - q.getY()) / (p.getX() - q.getX());

		return new MyDroite(-1 / a, milieu.y + (1 / a) * milieu.x);
	}
	
	public class MyDroite {

		double a, b;

		public MyDroite(double a, double b) {
			this.a = a;
			this.b = b;
		}

		public double getA() {
			return a;
		}

		public void setA(double a) {
			this.a = a;
		}

		public double getB() {
			return b;
		}

		public void setB(double b) {
			this.b = b;
		}
	}

}
