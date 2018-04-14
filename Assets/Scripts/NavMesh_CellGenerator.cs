using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.AI;
using UnityEngine.Rendering;

public enum InfluenceMode
{
	Propagation = 0,

}

public class Pair<T1, T2>
{
    public Pair(T1 f, T2 s)
    {
        first = f;
        second = s;
    }

    public T1 first { get; set; }
    public T2 second { get; set; }
}

public class NavMesh_CellGenerator : MonoBehaviour
{
    public GameObject navMeshVertPrefab;                        // Prefab for the object used to display NavMesh vertices
    public GameObject poissonDotPrefab;                         // Prefab for the object used to display Voronoi site locations
    public GameObject linePrefab;                               // Prefab for the LineRenderer used to display lines and edges
    public GameObject cellConnectionPrefab;                     // Prefab for the LineRenderer used to display edges between cells
    public List<Triangle> navMeshTris = new List<Triangle>();   // All Triangles on the NavMesh in the scene
    public float poissonTolerance = 0.5f;                       // What is the minimum distance apart each poisson point should be
    public float lineWidth = 0.01f;                             // How wide to draw all the lines in the scene
    public SortedCellList cells;                                // List of all existing cells within the scene -- will be used for the influence map calculations
	public float stepTime = 0.1f;
	float currentTime = 0.0f;
	public InfluenceMode mode = InfluenceMode.Propagation;
	public float DecayFactor = 0.2f;
	public float GrowthFactor = 0.3f;

	/// <summary>
	/// Represents one edge of a given Voronoi region
	/// </summary>
	public class Edge
    {
        public Vector3 start, end;      // Start and end of this edge
        public GameObject line;         // Line Renderer to visually display this edge in the scene

        public Edge(Vector3 _start, Vector3 _end, GameObject _linePrefab, float _lineWidth, Vector3 _normal, Color _color)
        {
            start = _start;
            end = _end;

            line = Instantiate(_linePrefab);
            LineRenderer lineRenderer = line.GetComponent<LineRenderer>();
            lineRenderer.SetPositions(new Vector3[2] { start, end });
            lineRenderer.enabled = true;
            lineRenderer.widthMultiplier = _lineWidth;
            lineRenderer.material.SetColor("_Color", _color);

            Transform colliderTransform = line.GetComponentInChildren<Transform>();

            colliderTransform.localPosition = ((end - start) / 2.0f) + start;
            colliderTransform.rotation = Quaternion.LookRotation(end - start, _normal);

            BoxCollider colliderBoxCollider = line.GetComponentInChildren<BoxCollider>();
            colliderBoxCollider.size = new Vector3(0.025f, 0.025f, (end - start).magnitude);
        }
    }

    /// <summary>
    /// Represents a single Voronoi region
    /// </summary>
    public class Cell
    {
        public List<Edge> edges;        // All edges belonging to this Voronoi region
        public Vector3 site;            // The origin of this Voronoi region

        public List<Pair<Cell, float>> adjacentCells; //A list of adjacent cells and the distances to them

		float newInfluenceValue;		// Newly calculated influence value
        float influenceValue;			// This Voronoi region's current influence (-1 to 1)
        GameObject dot;                 // Prefab PoissonDot object, should have a PoissonDot component attached
        Renderer dotRenderer;           // 
        GameObject linePrefab;          // 
		NavMesh_CellGenerator parent;	// Parent object, passed to the instantiated Dot GameObject

        private void ColorDot()
        {
            if (influenceValue < 0.0f)
                dotRenderer.material.SetColor("_Color", Math.Abs(influenceValue) * Color.blue);
            else
                dotRenderer.material.SetColor("_Color", influenceValue * Color.red);
        }

		private float GetInfluenceFromNeighbor(Pair<Cell, float> neighbor)
		{
			return neighbor.first.influenceValue * (float)(Math.Exp(-neighbor.second * parent.DecayFactor));
		}

		private float CalculateNewInfluence(float maxInfluence)
		{
			return (1.0f - parent.GrowthFactor) * influenceValue + parent.GrowthFactor * maxInfluence;
		}

        public Cell(Vector3 _site, Quaternion _rotation, GameObject _dotPrefab, GameObject _linePrefab, NavMesh_CellGenerator _parent)
        {
            site = _site;
            linePrefab = _linePrefab;
            influenceValue = 0.0f;

			parent = _parent;

			dot = Instantiate(_dotPrefab, site, _rotation);
            dotRenderer = dot.GetComponent<Renderer>();
            dotRenderer.material.SetColor("_Color", Color.white);

            adjacentCells = new List<Pair<Cell, float>>();
        }

		public void SetIndex(int index)
		{
			dot.GetComponent<PoissonDot>().SetParentAndIndex(parent, index);
		}

		public void SetInfluenceValue(float _influenceValue)
		{
			influenceValue = _influenceValue;
			ColorDot();
		}

		public void ApplyNewInfluence()
		{
			influenceValue = newInfluenceValue;
			ColorDot();
		}

		public void CalculateCurrentInfluence()
		{
			//List<float> influences = new List<float>(adjacentCells.Count);
			float maxInfluence = 0.0f;
			for(int i = 0; i < adjacentCells.Count; ++i)
			{
				float influence = GetInfluenceFromNeighbor(adjacentCells[i]);
				maxInfluence = Math.Abs(maxInfluence) < Math.Abs(influence) ? influence : maxInfluence;
				//influences[i] = GetInfluenceFromNeighbor(adjacentCells[i]);
			}

			newInfluenceValue = CalculateNewInfluence(maxInfluence);
		}
	}

    /// <summary>
    /// A single triangle from the NavMesh
    /// </summary>
    public class Triangle
    {
        public Vector3 v1, v2, v3;          // Vertices of this triangle
        public Vector3 normal;              // Normal of the plane spanning this triangle
        public GameObject l1, l2, l3;       // All Linerenderers associated with this triangle
        public float area;                  // Area of this triangle -- used in PoissonDiscDistribution
        public List<Cell> poissonCells;     // List of all Voronoi regions within this triangle
        public GameObject dotPrefab;        // Prefab for the object used to display Voronoi site locations
        public GameObject linePrefab;       // 
        public float lineWidth;             //

        /// <summary>
        /// Spawn a Line prefab and orient it's box collider properly
        /// </summary>
        /// <param name="line">Reference to the current line GameObject to initialize</param>
        /// <param name="start">Start point of the line</param>
        /// <param name="end">End point of the line</param>
        private void InitializeLine(ref GameObject line, Vector3 start, Vector3 end)
        {
            const float samplesPerLength = 2.0f;
            const float verticalOffset = 0.1f;

            line = Instantiate(linePrefab);
            LineRenderer lineRenderer = line.GetComponent<LineRenderer>();

            int samples = 2 + (int)(Vector3.Distance(start, end) * samplesPerLength);
            lineRenderer.positionCount = samples;
            Vector3[] linePositions = new Vector3[samples];
            for (int i = 0; i < samples; ++i)
            {
                linePositions[i] = Vector3.Lerp(start, end, (float)i / (samples - 1));
                NavMeshHit hit;
                NavMesh.SamplePosition(linePositions[i], out hit, 100, -1);
                linePositions[i] = hit.position + Vector3.up * verticalOffset;
            }
            lineRenderer.SetPositions(linePositions);

            //lineRenderer.startColor = Color.red;
            //lineRenderer.endColor = Color.red;
            lineRenderer.enabled = true;
            lineRenderer.widthMultiplier = lineWidth;

            Transform colliderTransform = line.GetComponentInChildren<Transform>();

            colliderTransform.localPosition = ((end - start) / 2.0f) + start;
            colliderTransform.rotation = Quaternion.LookRotation(end - start, normal);

            BoxCollider colliderBoxCollider = line.GetComponentInChildren<BoxCollider>();
            colliderBoxCollider.size = new Vector3(0.025f, 0.025f, (end - start).magnitude);
        }

        public Triangle(Vector3 _v1, Vector3 _v2, Vector3 _v3, GameObject _dotPrefab, GameObject _linePrefab, float _lineWidth)
        {
            v1 = _v1;
            v2 = _v2;
            v3 = _v3;

            dotPrefab = _dotPrefab;
            linePrefab = _linePrefab;
            lineWidth = _lineWidth;

            InitializeLine(ref l1, v1, v2);
            InitializeLine(ref l2, v2, v3);
            InitializeLine(ref l3, v3, v1);

            Vector3 a = v2 - v1;
            Vector3 b = v3 - v1;
            float aOnBProjectionLength = Vector3.Dot(a, b / b.magnitude);
            Vector3 aOnBProjection = b.normalized * aOnBProjectionLength;
            float height = (v2 - aOnBProjection).magnitude;
            float baseWidth = b.magnitude;
            area = (baseWidth * height) / 2.0f;

            normal = Vector3.Cross(a, b);
            normal.Normalize();
        }

        /// <summary>
        /// Returns a point somewhere within this triangle
        /// </summary>
        /// <returns></returns>
        public Vector3 RandomPoint()
        {
            float f1 = UnityEngine.Random.Range(0.0f, 1.0f);
            float f2 = UnityEngine.Random.Range(0.0f, 1.0f - f1);
            float f3 = 1.0f - f1 - f2;

            return (v1 * f1) + (v2 * f2) + (v3 * f3);
        }

        /// <summary>
        /// Assumes poissonCells are already set
        /// This method will create all Voronoi regions within this triangle
        /// </summary>
        public void CreateVoronoiRegions()
        {
            for (int outer = 0; outer < poissonCells.Count; ++outer)
            {
                for (int inner = 0; inner < poissonCells.Count; ++inner)
                {
                    if (inner == outer) continue;
                    Vector3 midPoint = (poissonCells[inner].site - poissonCells[outer].site) / 2.0f + poissonCells[outer].site;
                    Vector3 perp = Vector3.Cross(poissonCells[inner].site - poissonCells[outer].site, normal);

                    RaycastHit hit;
                    if (Physics.Raycast(midPoint, perp, out hit))
                    {
                        GameObject endPoint = Instantiate(dotPrefab, hit.point, new Quaternion());
                        endPoint.GetComponent<Renderer>().material.SetColor("_Color", Color.cyan);

                    }
                    if (Physics.Raycast(midPoint, -perp, out hit))
                    {
                        GameObject endPoint = Instantiate(dotPrefab, hit.point, new Quaternion());
                        endPoint.GetComponent<Renderer>().material.SetColor("_Color", Color.cyan);
                    }
                }
            }

            float oneEighth = 1.0f / 8.0f;
            for (int outer = 0; outer < poissonCells.Count; ++outer)
            {
                for (int inner = 0; inner < poissonCells.Count; ++inner)
                {
                    if (inner == outer) continue;
                    Vector3 midPoint = (poissonCells[inner].site - poissonCells[outer].site) / 2.0f + poissonCells[outer].site;
                    Vector3 perp = Vector3.Cross(poissonCells[inner].site - poissonCells[outer].site, normal);
                    perp = perp.normalized * oneEighth;
                    new Edge(midPoint, midPoint + perp, linePrefab, lineWidth, normal, Color.yellow);

                    new Edge(poissonCells[outer].site, midPoint, linePrefab, lineWidth, normal, inner > outer ? Color.green : Color.cyan);
                    new Edge(midPoint, midPoint + (normal / 4.0f), linePrefab, lineWidth, normal, Color.magenta);
                }
            }
        }
    }

    /// <summary>
    /// A list of Cells, sorted on the X axis
    /// </summary>
    public class SortedCellList
    {
        private List<Cell> cells = new List<Cell>();

        public void AddCell(Cell cell)
        {
            if (cells.Count == 0)
            {
                cells.Add(cell);
                return;
            }


            int i = cells.Count / 2;
            int offset = Mathf.Max(i / 2, 1);

            while (true)
            {
                if (i >= cells.Count)
                {
                    i = cells.Count;
                    if (cells[i - 1].site.x <= cell.site.x)
                    {
                        break;
                    }
                    --i;
                }
                else if (i <= 0)
                {
                    i = 0;
                    if (cells[0].site.x >= cell.site.x)
                    {
                        break;
                    }
                    ++i;
                }
                else if (cells[i].site.x >= cell.site.x && cells[i - 1].site.x <= cell.site.x)
                {
                    break;
                }
                if (cells[i - 1].site.x > cell.site.x)
                {
                    i -= offset;
                    offset = Mathf.Max(offset / 2, 1);
                }
                else
                {
                    i += offset;
                    offset = Mathf.Max(offset / 2, 1);
                }
            }
            cells.Insert(i, cell);

            //Debug info
            for (int j = 0; j < cells.Count - 1; ++j)
            {
                if (cells[j].site.x > cells[j + 1].site.x)
                {
                    Debug.Log("Not Sorted!!");
                }
            }
            //Debug.Log(cells.Count);
        }

        /// <summary>
        /// Finds the Cell with the X closest to value, without going below
        /// </summary>
        /// <param name="value"></param>
        /// <returns>The index of the element</returns>
        public int FindClosest(float value)
        {
            if (cells.Count == 0)
            {
                return 0;
            }
            int i = cells.Count / 2;
            int offset = Mathf.Max(i / 2, 1);

            while (true)
            {
                if (i >= cells.Count)
                {
                    i = cells.Count;
                    if (cells[i - 1].site.x <= value)
                    {
                        break;
                    }
                    --i;
                }
                else if (i <= 0)
                {
                    i = 0;
                    if (cells[0].site.x >= value)
                    {
                        break;
                    }
                    ++i;
                }
                else if (cells[i].site.x >= value && cells[i - 1].site.x <= value)
                {
                    break;
                }
                if (cells[i - 1].site.x >= value)
                {
                    i -= offset;
                    offset = Mathf.Max(offset / 2, 1);
                }
                else
                {
                    i += offset;
                    offset = Mathf.Max(offset / 2, 1);
                }
            }

            return i;
        }

        public Pair<int, int> GetRange(Vector3 center, float distance)
        {
            int low = FindClosest(center.x - distance);
            if (low != 0)
            {
                if (cells[low - 1].site.x >= center.x - distance)
                {
                    Debug.Log("Error! Point should have been included in the lower bound!");
                    Debug.Log(cells.Count);
                    Debug.Log(cells[low - 1].site.x + " vs. " + (center.x - distance));
                }
            }
            int high = FindClosest(center.x + distance);
            if (high != 0)
            {
                if (cells[high - 1].site.x > center.x + distance)
                {
                    //Debug.LogError("Error! Point should not have been included in the upper bound! (high too high)");
                    Debug.Log("Error! Point should not have been included in the upper bound! (high too high)");
                    Debug.Log(cells.Count);
                    Debug.Log(cells[high - 1].site.x + " vs. " + (center.x + distance));
                }
            }
            if (high != cells.Count)
            {
                if (cells[high].site.x <= center.x + distance)
                {
                    Debug.Log("Error! Point should not have been included in the upper bound! (high too low)");
                    Debug.Log(cells.Count);
                    Debug.Log(cells[high].site.x + " vs. " + (center.x + distance));
                }
            }

            return new Pair<int, int>(FindClosest(center.x - distance), FindClosest(center.x + distance));
        }

        public Cell GetCell(int index)
        {
            return cells[index];
        }

        public int GetCount()
        {
            return cells.Count;
        }

		/// <summary>
		/// For every cell contained within this list, calls cell.SetIndex(index) passing 
		/// the cell's position within the List as the index
		/// </summary>
		public void AssignIndices()
		{
			for (int i = 0; i < cells.Count - 1; ++i)
			{
				cells[i].SetIndex(i);
			}
		}

		/// <summary>
		/// Iterates through the cells and calls CalculateCurrentInfluence
		/// Once all cells have finished calculating influence, 
		/// it will iterate through the cells and call ApplyNewInfluence
		/// </summary>
		public void PropagateCells()
		{
			foreach(Cell cell in cells)
			{
				cell.CalculateCurrentInfluence();
			}

			foreach (Cell cell in cells)
			{
				cell.ApplyNewInfluence();
			}
		}
	}


	/// <summary>
	/// Given a Triangle and a tolerance, the triangle's cells member will be 
	/// populated with Cells that represent Voronoi regions
	/// </summary>
	/// <param name="triangle">Triangle to generate Poisson points within</param>
	/// <param name="tolerance">Minimum distance all Poisson points need to be away from one another</param>
	void PoissonDiscDistribution(Triangle triangle, float tolerance)
    {
        int numPoints = (int)(triangle.area * 100.0f) + 100;
        List<Vector3> points = new List<Vector3>();         // Initial points that are generated, enough to decently fill the triangle
        List<Vector3> pointsFinal = new List<Vector3>();    // Remaining points that have been filtered via the tolerance amount
        triangle.poissonCells = new List<Cell>();

        // Generate random points within the triangle
        for (int idx = 0; idx < numPoints; ++idx)
        {
            points.Add(triangle.RandomPoint());
            NavMeshHit hit;
            NavMesh.SamplePosition(points[points.Count - 1], out hit, 100, -1);
            points[points.Count - 1] = hit.position;

            //Here for stress testing
            //cells.AddCell(new Cell(points[points.Count - 1], transform.rotation, poissonDotPrefab, linePrefab));
            //cells.GetRange(points[points.Count - 1], poissonTolerance * 2);
        }

        // Test each point against the selected final points
        // If it meets the tolerance requirement, add it to pointsFinal
        // else, cull this point
        foreach (Vector3 point in points)
        {
            bool tooClose = false;
            Pair<int, int> bounds = cells.GetRange(point, poissonTolerance * 2);
            for (int i = bounds.first; i < bounds.second; ++i)
            {
                if (Vector3.Distance(point, cells.GetCell(i).site) < tolerance)
                {
                    tooClose = true;
                    break;
                }
            }

            if (!tooClose)
            {
                cells.AddCell(new Cell(point, transform.rotation, poissonDotPrefab, linePrefab, this));
                pointsFinal.Add(point);
            }
        }

        //// For all the selected points, generate Voronoi regions
        //foreach (Vector3 point in pointsFinal)
        //{
        //    triangle.poissonCells.Add(new Cell(point, transform.rotation, poissonDotPrefab, linePrefab));
        //}
        //
        //triangle.CreateVoronoiRegions();
    }

    bool PointWithinTriangleCircumference(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 v)
    {
        //Matrix4x4 matrix = Matrix4x4.zero;
        //matrix.m00 = A.x;
        //matrix.m01 = A.y;
        //matrix.m02 = A.x * A.x + A.y * A.y;
        //matrix.m03 = 1;
        //
        //matrix.m10 = B.x;
        //matrix.m11 = B.y;
        //matrix.m12 = B.x * B.x + B.y * B.y;
        //matrix.m13 = 1;
        //
        //matrix.m20 = C.x;
        //matrix.m21 = C.y;
        //matrix.m22 = C.x * C.x + C.y * C.y;
        //matrix.m23 = 1;
        //
        //matrix.m30 = D.x;
        //matrix.m31 = D.y;
        //matrix.m32 = D.x * D.x + D.y * D.y;
        //matrix.m33 = 1;

        ////
        //Algorithm for determining if a point is with the circumference of a triangle taken from: https://github.com/Bl4ckb0ne/delaunay-triangulation/blob/master/triangle.h
        float ab = (p1.x * p1.x) + (p1.z * p1.z);
        float cd = (p2.x * p2.x) + (p2.z * p2.z);
        float ef = (p3.x * p3.x) + (p3.z * p3.z);

        float circum_x = (ab * (p3.z - p2.z) + cd * (p1.z - p3.z) + ef * (p2.z - p1.z)) / (p1.x * (p3.z - p2.z) + p2.x * (p1.z - p3.z) + p3.x * (p2.z - p1.z)) / 2.0f;
        float circum_y = (ab * (p3.x - p2.x) + cd * (p1.x - p3.x) + ef * (p2.x - p1.x)) / (p1.z * (p3.x - p2.x) + p2.z * (p1.x - p3.x) + p3.z * (p2.x - p1.x)) / 2.0f;
        float circum_radius = Mathf.Sqrt(((p1.x - circum_x) * (p1.x - circum_x)) + ((p1.z - circum_y) * (p1.z - circum_y)));

        float dist = Mathf.Sqrt(((v.x - circum_x) * (v.x - circum_x)) + ((v.z - circum_y) * (v.z - circum_y)));
		////

		//D lies within the circumference of the circle
		return dist <= circum_radius;
    }

    void GetNeighbors(Cell cell)
    {
        Pair<int, int> bounds = cells.GetRange(cell.site, poissonTolerance * 10);

        if (bounds.first == bounds.second)
        {
            Debug.Log("No nearby points at all!!");
            return;
        }

        int nearest = -1;
        float nearestDist = Mathf.Infinity;
        List<Cell> possibleCells = new List<Cell>(bounds.second - bounds.first);
        List<Cell> allCells = new List<Cell>(bounds.second - bounds.first);
        for (int i = bounds.first; i < bounds.second; ++i)
        {
            float dist = Vector3.Distance(cell.site, cells.GetCell(i).site);
            if (dist < poissonTolerance || Mathf.Abs(cell.site.y - cells.GetCell(i).site.y) > poissonTolerance)
            {
                continue;
            }
            if (dist < poissonTolerance * 3)
            {
                possibleCells.Add(cells.GetCell(i));
                if (dist < nearestDist)
                {
                    nearest = possibleCells.Count - 1;
                    nearestDist = dist;
                }
            }
            allCells.Add(cells.GetCell(i));
        }

        if (possibleCells.Count == 0)
        {
            Debug.Log("No nearby points after distance check!!");
            return;
        }

        //Debug line to the nearest cell
        //GameObject l = Instantiate(linePrefab);
        //LineRenderer lr = l.GetComponent<LineRenderer>();
        //
        //int ss = 2 + (int)(Vector3.Distance(cell.site, possibleCells[nearest].site) * 3);
        //lr.positionCount = ss;
        //Vector3[] lps = new Vector3[ss];
        //for (int i = 0; i < ss; ++i)
        //{
        //	lps[i] = Vector3.Lerp(cell.site, possibleCells[nearest].site, (float)i / (ss - 1));
        //	NavMeshHit hit;
        //	NavMesh.SamplePosition(lps[i], out hit, 100, -1);
        //	lps[i] = hit.position + Vector3.up * 0.1f * ((float)i / (ss - 1));
        //}
        //lr.SetPositions(lps);
        //
        //lr.enabled = true;
        //lr.widthMultiplier = lineWidth;
        //lr.endWidth = lineWidth / 2;

        float nearestAngle = Vector3.SignedAngle(Vector3.right, possibleCells[nearest].site - cell.site, Vector3.up) + 180;

        //Bubble Sort of angles, low to high (measured clockwise from the nearest node)
        for (int i = 0; i < possibleCells.Count; ++i)
        {
            //Get the angle clockwise from the nearest node
            int lowestAngleIndex = i;
            float lowestAngle = nearestAngle - (Vector3.SignedAngle(Vector3.right, possibleCells[i].site - cell.site, Vector3.up) + 180);
            if (lowestAngle + 0.0001f < 0)
            {
                lowestAngle += 360;
            }
            //Find if there's a lower angle
            for (int j = i + 1; j < possibleCells.Count; ++j)
            {
                float angle = nearestAngle - (Vector3.SignedAngle(Vector3.right, possibleCells[j].site - cell.site, Vector3.up) + 180);
                if (angle + 0.0001f < 0)
                {
                    angle += 360;
                }

                if (angle < lowestAngle)
                {
                    lowestAngle = angle;
                    lowestAngleIndex = j;
                }
            }
            //Swap
            Cell temp = possibleCells[i];
            possibleCells[i] = possibleCells[lowestAngleIndex];
            possibleCells[lowestAngleIndex] = temp;
        }

        //Debug print
        //for (int i = 0; i < possibleCells.Count; ++i)
        //{
        //	float angle = nearestAngle - (Vector3.SignedAngle(Vector3.right, possibleCells[i].site - cell.site, Vector3.up) + 180);
        //	if (angle < 0)
        //	{
        //		angle += 360;
        //	}
        //	Debug.Log(i + " index has angle " + angle);
        //}
        //if (Vector3.Distance(possibleCells[0].site, cell.site) != nearestDist)
        //{
        //	Debug.Log(Vector3.Distance(possibleCells[0].site, cell.site) + " vs. nearest of " + nearestDist);
        //}

        List<Cell> finalCells = new List<Cell>();
        finalCells.Add(possibleCells[0]);
        for (int i = 1; i < possibleCells.Count; ++i)
        {
            bool pointIsNeighbor = true;
            for (int j = 0; j < allCells.Count; ++j)
            {
                if (finalCells[finalCells.Count - 1] == allCells[j] || possibleCells[i] == allCells[j])
                {
                    continue;
                }
                if (PointWithinTriangleCircumference(finalCells[finalCells.Count - 1].site, cell.site, possibleCells[i].site, allCells[j].site))
                {
                    pointIsNeighbor = false;
                    break;
                }
            }
            if (pointIsNeighbor)
            {
                finalCells.Add(possibleCells[i]);
            }
        }


        foreach (Cell c in finalCells)
        {
            const float samplesPerLength = 2.0f;
            const float verticalOffset = 0.0f;

            GameObject line = Instantiate(cellConnectionPrefab);
            LineRenderer lineRenderer = line.GetComponent<LineRenderer>();

            int samples = 2 + (int)(Vector3.Distance(cell.site, c.site) * samplesPerLength);
            lineRenderer.positionCount = samples;
            Vector3[] linePositions = new Vector3[samples];
            for (int i = 0; i < samples; ++i)
            {
                linePositions[i] = Vector3.Lerp(cell.site, c.site, (float)i / (samples - 1));
                NavMeshHit hit;
                NavMesh.SamplePosition(linePositions[i], out hit, 100, -1);
                //float angle = nearestAngle - (Vector3.SignedAngle(Vector3.right, c.site - cell.site, Vector3.up) + 180);
                linePositions[i] = hit.position + Vector3.up * verticalOffset * ((float)i / (samples - 1));
            }
            lineRenderer.SetPositions(linePositions);

            lineRenderer.enabled = true;
            lineRenderer.widthMultiplier = lineWidth;
            lineRenderer.endWidth = lineWidth / 4;
            cell.adjacentCells.Add(new Pair<Cell, float>(c, Vector3.Distance(c.site, cell.site))); //Add
        }
    }

	/// <summary>
	/// Given an index into cells, sets that Cell's influence value
	/// </summary>
	/// <param name="index">Index into cells List</param>
	/// <param name="influence">Influence value to set this cell to</param>
	public void PropagationFromCell(int index, float influence)
	{
		cells.GetCell(index).SetInfluenceValue(influence);
	}

    // Use this for initialization
    void Start()
    {
        cells = new SortedCellList();

        // Get all NavMesh information
        NavMeshTriangulation tris = NavMesh.CalculateTriangulation();

        // Spawn vertex prefabs at the vertices of the NavMesh
        foreach (Vector3 vert in tris.vertices)
        {
            Instantiate(navMeshVertPrefab, vert, transform.rotation);
        }

        // Using the indicies of the NavMesh, create Triangles and call PoissonDiscDistribution with it to create
        // Voronoi regions
        for (int index = 0; index < tris.indices.Length; index += 3)
        {
            Triangle tri = new Triangle(
                tris.vertices[tris.indices[index + 0]],
                tris.vertices[tris.indices[index + 1]],
                tris.vertices[tris.indices[index + 2]],
                poissonDotPrefab, linePrefab, lineWidth
            );
            navMeshTris.Add(tri);
            PoissonDiscDistribution(tri, poissonTolerance);
        }

        for (int i = 0; i < cells.GetCount(); ++i)
        {
            GetNeighbors(cells.GetCell(i));
        }

		cells.AssignIndices();
    }

    // Update is called once per frame
    void Update()
    {
		if (mode.Equals(InfluenceMode.Propagation))
		{
			currentTime += Time.deltaTime;

			if (currentTime >= stepTime)
			{
				cells.PropagateCells();
			}
		}
    }
}
