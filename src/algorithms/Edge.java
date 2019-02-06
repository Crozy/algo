package algorithms;

import java.awt.Point;

public class Edge {
	private Point a;
	private Point b;
	private double poids;
	
	public Edge(Point a, Point b) {
		this.a = a;
		this.b = b;
		this.poids = DefaultTeam.distanceSQRT(a, b);
	}
	public Point getA() {
		return a;
	}
	public void setA(Point a) {
		this.a = a;
	}
	public Point getB() {
		return b;
	}
	public void setB(Point b) {
		this.b = b;
	}
	public double getPoids() {
		return poids;
	}
	public void setPoids(double poids) {
		this.poids = poids;
	}
		
	@Override
	public String toString() {
		return poids + " x: " + a.x + " y: " + a.y + " x: " + b.x + " y: " + b.y +"\n";
	}


}
