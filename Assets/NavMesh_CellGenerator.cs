using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.AI;
using UnityEngine.Rendering;


public class NavMesh_CellGenerator : MonoBehaviour {

    public GameObject navMeshVertPrefab;                        // Prefab for the object used to display NavMesh vertices
    public GameObject navMeshPoissonDotPrefab;                  // Prefab for the object used to display Voronoi site locations
    public GameObject lineRendererPrefab;                       // Prefab for the LineRenderer used to display lines and edges
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
        public LineRenderer line;       // Line Renderer to visually display this edge in the scene

        public Edge(Vector3 _start, Vector3 _end, GameObject linePrefab, float _lineWidth)
        {
            start = _start;
            end = _end;

            line = Instantiate(linePrefab).GetComponent<LineRenderer>();
            line.SetPositions(new Vector3[2] { start, end });
            line.enabled = true;
            line.widthMultiplier = _lineWidth;
        }
    }

    /// <summary>
    /// Represents a single Voronoi region
    /// </summary>
    public class Cell
    {
        public List<Edge> edges;        // All edges belonging to this Voronoi region
        public Vector3 site;            // The origin of this Voronoi region
        public GameObject linePrefab;   // Line prefab to send to Edge objects
        public float lineWidth;         // How wide to render this cell's lines

        public float influenceValue;    // This Voronoi region's current influence (-1 to 1)
        public GameObject dot;          // Dot prefab 
        Renderer dotRenderer;
        

        public Cell(GameObject _linePrefab, float _lineWidth, GameObject dotPrefab, Vector3 _site, Quaternion rotation)
        {
            linePrefab = _linePrefab;
            lineWidth = _lineWidth;
            site = _site;

            dot = Instantiate(dotPrefab, site, rotation);
            dotRenderer = dot.GetComponent<Renderer>();
            dotRenderer.material.SetColor("_Color", Color.blue);
        }

        void ColorDot()
        {
            if (influenceValue < 0.0f)
                dotRenderer.material.SetColor("_Color", influenceValue * Color.blue);
            else
                dotRenderer.material.SetColor("_Color", influenceValue * Color.red);
        }
    }

    /// <summary>
    /// A single triangle from the NavMesh
    /// </summary>
    public class Triangle
    {
        public Vector3 v1, v2, v3;          // Vertices of this triangle
        public LineRenderer l1, l2, l3;     // All Linerenderers associated with this triangle
        public float lineWidth;             // How wide to draw all lines of this triangle
        public float area;                  // Area of this triangle -- used in PoissonDiscDistribution
        public List<Cell> poissonCells;     // List of all Voronoi regions within this triangle

        public Triangle(Vector3 _v1, Vector3 _v2, Vector3 _v3, GameObject linePrefab, float _lineWidth)
        {
            v1 = _v1;
            v2 = _v2;
            v3 = _v3;

            lineWidth = _lineWidth;

            l1 = Instantiate(linePrefab).GetComponent<LineRenderer>();
            l1.SetPositions(new Vector3[2] { v1, v2 });
            l1.startColor = Color.red;
            l1.endColor = Color.red;
            l1.enabled = true;
            l1.widthMultiplier = lineWidth;

            l2 = Instantiate(linePrefab).GetComponent<LineRenderer>();
            l2.SetPositions(new Vector3[2] { v2, v3 });
            l2.startColor = Color.red;
            l2.endColor = Color.red;
            l2.enabled = true;
            l2.widthMultiplier = lineWidth;

            l3 = Instantiate(linePrefab).GetComponent<LineRenderer>();
            l3.SetPositions(new Vector3[2] { v3, v1 });
            l3.startColor = Color.red;
            l3.endColor = Color.red;
            l3.enabled = true;
            l3.widthMultiplier = lineWidth;

            Vector3 a = v2 - v1;
            Vector3 b = v3 - v1;
            float aOnBProjectionLength = Vector3.Dot(a, b / b.magnitude);
            Vector3 aOnBProjection = b.normalized * aOnBProjectionLength;
            float height = (v2 - aOnBProjection).magnitude;
            float baseWidth = b.magnitude;
            area = (baseWidth * height) / 2.0f;
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
            triangle.poissonCells.Add(new Cell(lineRendererPrefab, lineWidth, navMeshPoissonDotPrefab, point, transform.rotation));
        }
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
                lineRendererPrefab,
                lineWidth
            );
            navMeshTris.Add(tri);
            PoissonDiscDistribution(tri, poissonTolerance);
        }
    }
	
	// Update is called once per frame
	void Update () {

    }
}
