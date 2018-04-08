using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.AI;
using UnityEngine.Rendering;


public class NavMesh_CellGenerator : MonoBehaviour {

    public GameObject navMeshVertPrefab;                        // Prefab for the object used to display NavMesh vertices
    public GameObject poissonDotPrefab;                         // Prefab for the object used to display Voronoi site locations
    public GameObject linePrefab;                               // Prefab for the LineRenderer used to display lines and edges
    public List<Triangle> navMeshTris = new List<Triangle>();   // All Triangles on the NavMesh in the scene
    public float poissonTolerance = 0.5f;                       // What is the minimum distance apart each poisson point should be
    public float lineWidth = 0.01f;                             // How wide to draw all the lines in the scene
    public List<Cell> cells;                                    // List of all existing cells within the scene -- will be used for the influence map calculations

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

        public float influenceValue;    // This Voronoi region's current influence (-1 to 1)
        GameObject dot;                 // 
        Renderer dotRenderer;           // 
        GameObject linePrefab;          // 
        
        private void ColorDot()
        {
            if (influenceValue < 0.0f)
                dotRenderer.material.SetColor("_Color", influenceValue * Color.blue);
            else
                dotRenderer.material.SetColor("_Color", influenceValue * Color.red);
        }

        public Cell(Vector3 _site, Quaternion _rotation, GameObject _dotPrefab, GameObject _linePrefab)
        {
            site = _site;
            linePrefab = _linePrefab;
            influenceValue = 0.0f;

            dot = Instantiate(_dotPrefab, site, _rotation);
            dotRenderer = dot.GetComponent<Renderer>();
            dotRenderer.material.SetColor("_Color", Color.white);
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
            line = Instantiate(linePrefab);
            LineRenderer lineRenderer = line.GetComponent<LineRenderer>();
            lineRenderer.SetPositions(new Vector3[2] { start, end });
            //lineRenderer.startColor = Color.red;
            //lineRenderer.endColor = Color.red;
            lineRenderer.enabled = true;
            lineRenderer.widthMultiplier = lineWidth;

            Transform colliderTransform = line.GetComponentInChildren<Transform>();
            
            colliderTransform.localPosition = ((end - start) / 2.0f) + start;
            colliderTransform.rotation = Quaternion.LookRotation(end-start, normal);

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
            float f1 = Random.Range(0.0f, 1.0f);
            float f2 = Random.Range(0.0f, 1.0f - f1);
            float f3 = 1.0f - f1 - f2;

            return (v1 * f1) + (v2 * f2) + (v3 * f3);
        }

        /// <summary>
        /// Assumes poissonCells are already set
        /// This method will create all Voronoi regions within this triangle
        /// </summary>
        public void CreateVoronoiRegions()
        {
            for(int outer = 0; outer< poissonCells.Count; ++outer)
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
    /// Given a Triangle and a tolerance, the triangle's cells member will be 
    /// populated with Cells that represent Voronoi regions
    /// </summary>
    /// <param name="triangle">Triangle to generate Poisson points within</param>
    /// <param name="tolerance">Minimum distance all Poisson points need to be away from one another</param>
    void PoissonDiscDistribution(Triangle triangle, float tolerance)
    {
        int numPoints = (int)(triangle.area * 25.0f) + 100;
        List<Vector3> points = new List<Vector3>();         // Initial points that are generated, enough to decently fill the triangle
        List<Vector3> pointsFinal = new List<Vector3>();    // Remaining points that have been filtered via the tolerance amount
        triangle.poissonCells = new List<Cell>();

        // Generate random points within the triangle
        for (int idx = 0; idx < numPoints; ++idx)
        {
            points.Add(triangle.RandomPoint());
        }

        // Test each point against the selected final points
        // If it meets the tolerance requirement, add it to pointsFinal
        // else, cull this point
        foreach (Vector3 point in points)
        {
            bool tooClose = false;
            foreach (Vector3 finalpoint in pointsFinal)
            {
                if (Vector3.Distance(point, finalpoint) < tolerance)
                {
                    tooClose = true;
                    break;
                }
            }

            if (!tooClose)
                pointsFinal.Add(point);
        }

        // For all the selected points, generate Voronoi regions
        foreach (Vector3 point in pointsFinal)
        {
            triangle.poissonCells.Add(new Cell(point, transform.rotation, poissonDotPrefab, linePrefab));
        }

        triangle.CreateVoronoiRegions();
    }

    // Use this for initialization
    void Start () {
        cells = new List<Cell>();

        // Get all NavMesh information
        NavMeshTriangulation tris = NavMesh.CalculateTriangulation();

        // Spawn vertex prefabs at the vertices of the NavMesh
        foreach (Vector3 vert in tris.vertices)
        {
            Instantiate(navMeshVertPrefab, vert, transform.rotation);
        }

        // Using the indicies of the NavMesh, create Triangles and call PoissonDiscDistribution with it to create
        // Voronoi regions
        for(int index = 0; index < tris.indices.Length; index+=3)
        {
            Triangle tri = new Triangle(
                tris.vertices[tris.indices[index+0]],
                tris.vertices[tris.indices[index+1]],
                tris.vertices[tris.indices[index+2]],
                poissonDotPrefab, linePrefab, lineWidth
            );
            navMeshTris.Add(tri);
            PoissonDiscDistribution(tri, poissonTolerance);
        }
    }
	
	// Update is called once per frame
	void Update () {
    }
}
