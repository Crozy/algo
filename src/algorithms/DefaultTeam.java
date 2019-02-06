package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import algorithms.Tools.MyDroite;

public class DefaultTeam {
	double[][] D;
	int[][] path;
	public static final double MYTHICAL_NUMBER = 1664;

	public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold,
			ArrayList<Point> hitPoints) {
		double oldscore, newscore;
		path = calculShortestPaths(points, edgeThreshold);

		ArrayList<Edge> result = tme5(points, edgeThreshold, hitPoints);
		oldscore = getScore(result);
		result = tme4Barycentre(getPointsFromArretes(result), points,
				hitPoints, edgeThreshold); //
		newscore = getScore(result);

		while (newscore < oldscore) {
			oldscore = newscore;
			result = tme4Barycentre(getPointsFromArretes(result), points,
					hitPoints, edgeThreshold);
			newscore = getScore(result);
		}
		return toTree2D(result.get(0).getA(), result);
	}
	
	public ArrayList<Point> calculSteinerPoints(ArrayList<Point> points, int edgeThreshold,
			ArrayList<Point> hitPoints) {
		double oldscore, newscore;
		path = calculShortestPaths(points, edgeThreshold);

		ArrayList<Edge> result = tme5(points, edgeThreshold, hitPoints);
		oldscore = getScore(result);
		result = tme4Barycentre(getPointsFromArretes(result), points,
				hitPoints, edgeThreshold); //
		newscore = getScore(result);
		
		
		while (newscore < oldscore) {
			oldscore = newscore;
			result = tme4Barycentre(getPointsFromArretes(result), points,
					hitPoints, edgeThreshold);
			newscore = getScore(result);
		}
		
		Set<Point> res = new HashSet<Point>();
		for (Edge edge : result) {
			res.add(edge.getA());
			res.add(edge.getB());
		}
		
		return new ArrayList<Point>(res);
	}

	public ArrayList<Edge> getKruskal(ArrayList<Point> points) {
		ArrayList<Edge> listTri = allArrete(points);
		ArrayList<Edge> res = new ArrayList<>();
		int i = 0;

		while (i < listTri.size()) {
			if (!detectCycle(res, listTri.get(i)))
				res.add(listTri.get(i));
			i++;
		}
		return res;
	}

	public Tree2D calculSteinerBudget(ArrayList<Point> points,
			int edgeThreshold, ArrayList<Point> hitPoints) {

		path = calculShortestPaths(points, edgeThreshold);
		ArrayList<Point> res = new ArrayList<Point>();
		res.add(hitPoints.remove(0));
		ArrayList<Point> hitMinus = new ArrayList<Point>(hitPoints);

		boolean find = false;
		double budget = 0;
		SpecialPoint spe = null;

		for (int i = 0; i < hitPoints.size(); i++) {
			while (budget <= MYTHICAL_NUMBER) {
				spe = pointLePlusProche(res, hitMinus, points);
				if (budget + spe.dist <= MYTHICAL_NUMBER) {
					budget += spe.dist;
					res.add(spe.p);
					hitMinus.remove(spe.p);
				} else if(budget+spe.dist > MYTHICAL_NUMBER) {
					find = true;
					break;
				}
			}

			if (find)
				break;
		}

		return toTree2D(res.get(0), tme5(points, edgeThreshold, res));
	}

	public int[][] calculShortestPaths(ArrayList<Point> points,
			int edgeThreshold) {
		int[][] paths = new int[points.size()][points.size()];
		D = new double[points.size()][points.size()];
		double dist;

		for (int i = 0; i < points.size(); i++) {
			for (int j = 0; j < points.size(); j++) {
				dist = distanceSQRT(points.get(i), points.get(j));
				if (dist <= edgeThreshold)
					D[i][j] = 1;
				else if (i == j)
					D[i][j] = 0;
				else
					D[i][j] = Double.POSITIVE_INFINITY;
				paths[i][j] = j;
			}
		}

		for (int k = 1; k < points.size(); k++) {
			for (int i = 0; i < points.size(); i++) {
				for (int j = 0; j < points.size(); j++) {
					if (D[i][k] + D[k][j] < D[i][j]) {
						D[i][j] = D[i][k] + D[k][j];
						paths[i][j] = paths[i][k];
					}
				}
			}
		}

		return paths;
	}

	public ArrayList<Edge> tme5(ArrayList<Point> points, int edgeThreshold,
			ArrayList<Point> hitPoints) {
		ArrayList<Edge> graphecomplet = allEdges(hitPoints, points);
		ArrayList<Edge> res = new ArrayList<Edge>();
		Edge tmp;
		int i = 0;

		while (i < graphecomplet.size()) {
			if (!detectCycle(res, graphecomplet.get(i)))
				res.add(graphecomplet.get(i));
			i++;
		}
		int indexa = 0;
		int indexb = 0;
		ArrayList<Edge> H = new ArrayList<Edge>();
		for (Edge edge : res) {
			Point a = edge.getA();
			Point b = edge.getB();
			indexb = points.indexOf(b);
			indexa = points.indexOf(a);
			while (indexa != path[indexa][indexb]) {
				tmp = new Edge(a, points.get(path[indexa][indexb]));
				H.add(tmp);
				a = points.get(path[indexa][indexb]);
				indexa = points.indexOf(a);
			}
		}

		return H;
	}

	public ArrayList<Edge> tme4Barycentre(ArrayList<Point> points,
			ArrayList<Point> bluePoints, ArrayList<Point> hitpoints,
			int edgeTreshold) {
		ArrayList<Edge> listTri = allArrete(points);
		ArrayList<Edge> res = new ArrayList<Edge>();
		Point barycentre, a, b, c;
		barycentre = a = b = c = null;
		Edge a1, a2;
		double oldscore, newscore, poidsar, poidsbar;
		int i = 0;
		while (i < listTri.size()) {
			if (!detectCycle(res, listTri.get(i)))
				res.add(listTri.get(i));
			i++;
		}
		oldscore = getScore(res);

		for (i = 0; i < res.size(); i++) {
			for (int j = 0; j < res.size(); j++) {
				a1 = res.get(i);
				a2 = res.get(j);
				if (a1.equals(a2))
					continue;

				if (a1.getA().equals(a2.getA())) {
					barycentre = getBarycentre(a1.getB(), a1.getA(), a2.getB());
					a = a1.getA();
					b = a1.getB();
					c = a2.getB();
				} else if (a1.getA().equals(a2.getB())) {
					barycentre = getBarycentre(a1.getB(), a1.getA(), a2.getA());
					b = a1.getB();
					a = a1.getA();
					c = a2.getA();
				} else if (a1.getB().equals(a2.getA())) {
					barycentre = getBarycentre(a1.getA(), a1.getB(), a2.getB());
					b = a1.getA();
					a = a1.getB();
					c = a2.getB();
				} else if (a1.getB().equals(a2.getB())) {
					barycentre = getBarycentre(a1.getA(), a1.getB(), a2.getA());
					b = a1.getA();
					a = a1.getB();
					c = a2.getA();
				}
				if (barycentre == null)
					continue;

				barycentre = pointLePlusProche(barycentre, bluePoints);
				poidsar = D[bluePoints.indexOf(a1.getA())][bluePoints
						.indexOf(a1.getB())]
						+ D[bluePoints.indexOf(a2.getA())][bluePoints
								.indexOf(a2.getB())];
				int indexbar = bluePoints.indexOf(barycentre);
				poidsbar = D[indexbar][bluePoints.indexOf(a)]
						+ D[indexbar][bluePoints.indexOf(b)]
						+ D[indexbar][bluePoints.indexOf(c)];
				newscore = getScore(res) + poidsbar - poidsar;

				if (newscore < oldscore) {
					hitpoints.add(barycentre);
					res = tme5(bluePoints, edgeTreshold, hitpoints);
					oldscore = newscore;
					newscore = getScore(res);
				}
				barycentre = null;
			}
		}
		return res;
	}

	public ArrayList<Edge> tme4FermatToricelli(ArrayList<Point> points,
			ArrayList<Point> bluePoints, ArrayList<Point> hitpoints,
			int edgeThreshold) {
		ArrayList<Edge> listTri = allArrete(points);
		ArrayList<Edge> res = new ArrayList<Edge>();
		Point fermatToricelli, a, b, c;
		fermatToricelli = a = b = c = null;
		Edge a1, a2;
		double oldscore, newscore, poidsar, poidsbar;
		int i = 0;

		while (res.size() < (points.size() - 1)) {
			if (!detectCycle(res, listTri.get(i)))
				res.add(listTri.get(i));
			i++;
		}
		oldscore = getScore(res);

		for (i = 0; i < res.size(); i++) {
			for (int j = 0; j < res.size(); j++) {
				a1 = res.get(i);
				a2 = res.get(j);
				if (a1.equals(a2))
					continue;

				if (a1.getA().equals(a2.getA())) {
					fermatToricelli = getFermatToricelli(a1.getB(), a1.getA(),
							a2.getB());
					a = a1.getA();
					b = a1.getB();
					c = a2.getB();
				} else if (a1.getA().equals(a2.getB())) {
					fermatToricelli = getFermatToricelli(a1.getB(), a1.getA(),
							a2.getA());
					b = a1.getB();
					a = a1.getA();
					c = a2.getA();
				} else if (a1.getB().equals(a2.getA())) {
					fermatToricelli = getFermatToricelli(a1.getA(), a1.getB(),
							a2.getB());
					b = a1.getA();
					a = a1.getB();
					c = a2.getB();
				} else if (a1.getB().equals(a2.getB())) {
					fermatToricelli = getFermatToricelli(a1.getA(), a1.getB(),
							a2.getA());
					b = a1.getA();
					a = a1.getB();
					c = a2.getA();
				}

				if (fermatToricelli == null)
					continue;

				fermatToricelli = pointLePlusProche(fermatToricelli, bluePoints);
				poidsar = D[bluePoints.indexOf(a1.getA())][bluePoints
						.indexOf(a1.getB())]
						+ D[bluePoints.indexOf(a2.getA())][bluePoints
								.indexOf(a2.getB())];
				int indexbar = bluePoints.indexOf(fermatToricelli);
				poidsbar = D[indexbar][bluePoints.indexOf(a)]
						+ D[indexbar][bluePoints.indexOf(b)]
						+ D[indexbar][bluePoints.indexOf(c)];
				newscore = getScore(res) + poidsbar - poidsar;

				if (newscore < oldscore) {
					hitpoints.add(fermatToricelli);
					res = tme5(bluePoints, edgeThreshold, hitpoints);
					oldscore = newscore;
					newscore = getScore(res);
				}
				fermatToricelli = null;
			}
		}

		return res;
	}

	public Point getFermatToricelli(Point a, Point b, Point c) {
		int x, y;

		Tools tools = new Tools();
		MyDroite d1 = tools.droitePerpendiculairePassantParMilieu(a, b);
		MyDroite d2 = tools.droitePerpendiculairePassantParMilieu(b, c);

		x = (int) ((d2.getB() - d1.getB()) / (d1.getA() - d2.getA()));
		y = (int) (d1.getA() * x + d1.getB());

		return new Point(x, y);
	}

	public static double distanceSQRT(Point a, Point b) {
		return Math.sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
	}

	public Tree2D toTree2D(Point res, ArrayList<Edge> list) {
		if (list.isEmpty())
			return new Tree2D(res, new ArrayList<Tree2D>());
		ArrayList<Tree2D> fils = new ArrayList<Tree2D>();
		ArrayList<Edge> edgeFromRes = new ArrayList<Edge>();
		for (Edge e : list) {
			if (e.getA().getX() == res.getX() && e.getA().getY() == res.getY()
					|| e.getB().getX() == res.getX()
					&& e.getB().getY() == res.getY()) {
				edgeFromRes.add(e);
			}
		}
		list.removeAll(edgeFromRes);
		for (Edge e : edgeFromRes) {
			if (e.getA().getX() == res.getX() && e.getA().getY() == res.getY()) {
				fils.add(toTree2D(e.getB(), list));
			} else if (e.getB().getX() == res.getX()
					&& e.getB().getY() == res.getY()) {
				fils.add(toTree2D(e.getA(), list));
			}
		}
		return new Tree2D(res, fils);
	}

	public ArrayList<Edge> allEdges(ArrayList<Point> points,
			ArrayList<Point> indexlist) {
		ArrayList<Edge> edges = new ArrayList<Edge>();

		for (int i = 0; i < points.size(); i++) {
			for (int j = i + 1; j < points.size(); j++) {
				int indexi = indexlist.indexOf(points.get(i));
				int indexj = indexlist.indexOf(points.get(j));
				Edge toadd = new Edge(points.get(i), points.get(j));
				toadd.setPoids(D[indexi][indexj]);
				edges.add(toadd);
			}
		}
		edges.sort(new EdgeComparator());

		return edges;
	}

	class EdgeComparator implements Comparator<Edge> {

		@Override
		public int compare(Edge o1, Edge o2) {
			if (o1.getPoids() < o2.getPoids())
				return -1;
			else if (o1.getPoids() > o2.getPoids())
				return 1;
			else
				return 0;
		}

	}

	public ArrayList<Edge> allArrete(ArrayList<Point> points) {
		ArrayList<Edge> edges = new ArrayList<Edge>();

		for (int i = 0; i < points.size(); i++)
			for (int j = i + 1; j < points.size(); j++)
				edges.add(new Edge(points.get(i), points.get(j)));
		edges.sort(new EdgeComparator());

		return edges;
	}

	public boolean detectCycle(ArrayList<Edge> list, Edge toAdd) {
		int i = 1;
		ArrayList<Edge> listCpy = new ArrayList<Edge>();
		listCpy.addAll(list);
		listCpy.add(toAdd);
		HashMap<Point, Integer> etiquette = new HashMap<Point, Integer>();
		for (Edge edge : listCpy) { // Donner une etiquette a tout les
										// points
			etiquette.put(edge.getA(), i++);
			etiquette.put(edge.getB(), i++);
		}
		for (Edge edge : listCpy) {
			int etiqA = etiquette.get(edge.getA());
			int etiqB = etiquette.get(edge.getB());
			if (etiqA == etiqB) {
				return true;
			}
			for (Point x : etiquette.keySet()) {
				if (etiquette.get(x) == etiqA) {
					etiquette.put(x, etiqB);
				}
			}
		}
		return false;
	}

	public ArrayList<Point> getPointsFromArretes(ArrayList<Edge> edges) {
		Set<Point> points = new HashSet<Point>();
		for (Edge edge : edges) {
			points.add(edge.getA());
			points.add(edge.getB());
		}
		return new ArrayList<Point>(points);
	}

	public double getScore(ArrayList<Edge> ar) {
		double res = 0.0;
		for (Edge a : ar)
			res += a.getPoids();
		return res;
	}

	public Point getBarycentre(Point a, Point b, Point c) {
		int xc = (a.x + b.x + c.x) / 3;
		int yc = (a.y + b.y + c.y) / 3;
		return new Point(xc, yc);
	}

	public Point getBarycentre(Point a, Point b, Point c, Point d) {
		int xc = (a.x + b.x + c.x + d.x) / 4;
		int yc = (a.y + b.y + c.y + d.y) / 4;
		return new Point(xc, yc);
	}

	private Point pointLePlusProche(Point barycentre,
			ArrayList<Point> bluePoints) {
		double dist = Double.MAX_VALUE, tmp;
		Point graal = null;

		for (Point point : bluePoints) {
			tmp = distanceSQRT(barycentre, point);
			if (tmp > 0 && tmp < dist) {
				graal = point;
				dist = tmp;
			}
		}
		return graal;
	}


	public SpecialPoint pointLePlusProche(ArrayList<Point> res,
			ArrayList<Point> hitMinus, ArrayList<Point> points) {
		Point proche = null;
		double minDist = Double.MAX_VALUE;

		for (Point p : res) {
			for (Point q : hitMinus) {
				if (D[points.indexOf(p)][points.indexOf(q)] < minDist) {
					minDist = D[points.indexOf(p)][points.indexOf(q)];
					proche = q;
				}
			}
		}
		return new SpecialPoint(proche, minDist);
	}

	/*
	 * Holly JAVA, i can't return a tuple.
	 */
	class SpecialPoint {
		Point p;
		double dist;

		public SpecialPoint(Point p, double dist) {
			this.p = p;
			this.dist = dist;
		}
	}
}