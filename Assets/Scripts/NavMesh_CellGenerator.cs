using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.AI;
using UnityEngine.Rendering;

public enum InfluenceMode
{
	Propagation = 0,			// Propagation as normal, click a cell to "drop" influence onto the map
	NormalizedPropagation,		// Same as above but values are normalized with max influence values
	FieldOfView,				// Allows the "player" to paint the cells via the line of sight box (Field of View)
	VisibleToSpot,				// How many other cells can "see" this one
	OpennessClosestWall,        // How close is this cell to a wall

	NumValues
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

/// <summary>
/// Represents a single Voronoi region
/// </summary>
public class Cell
{
	public Vector3 site;            // The origin of this Voronoi region

	public List<Pair<Cell, float>> adjacentCells; //A list of adjacent cells and the distances to them

	float currentFadeTime;          // Current time from MaxFadeTimer to 0 used to interpolate influence color
	float newInfluenceValue;        // Newly calculated influence value
	float influenceValue;           // This Voronoi region's current influence (-1 to 1)
	float wallInfluence;			// Stored influence value for closeness to a wall
	GameObject dot;                 // PoissonDot object, should have a PoissonDot component attached
	Renderer dotRenderer;           // Renderer component from the PoissonDot object spawned upon this cell's creation
	NavMesh_CellGenerator parent;   // Parent object, passed to the instantiated Dot GameObject

	/// <summary>
	/// Given a t value, performs linear interpolation between the color white and the given color
	/// </summary>
	/// <param name="t">Value between 0 and 1</param>
	/// <param name="otherColor">Color to interpolate with</param>
	/// <returns></returns>
	private Color LerpColor(float t, Color otherColor)
	{
		return (1.0f - t) * Color.white + t * otherColor;
	}

	/// <summary>
	/// If this cell's influenceValue is greater than or equal to 0, calls LerpColor with the 
	/// influenceValue and Red
	/// Else, calls LerpColor with the absolute influenceValue and Blue
	/// 
	/// Either result in this cell's dotRenderer having it's material color changed to the result of 
	/// LerpColor
	/// </summary>
	private void ColorDot()
	{
		if (influenceValue < 0.0f)
			dotRenderer.material.SetColor("_Color", LerpColor(Math.Abs(influenceValue), Color.blue));
		else
			dotRenderer.material.SetColor("_Color", LerpColor(influenceValue, Color.red));
	}

	/// <summary>
	/// For a given neighbor cell, calculates it's influence value using a decay formula and the 
	/// NavMesh_CellGenerator's DecayFactor
	/// </summary>
	/// <param name="neighbor">Pair containing the neighbor Cell and the distance from that cell to this one</param>
	/// <returns>The newly intrpolated influence value from the neighbor cell</returns>
	private float GetInfluenceFromNeighbor(Pair<Cell, float> neighbor)
	{
		return neighbor.first.influenceValue * (float)(Math.Exp(-neighbor.second * parent.DecayFactor));
	}

	/// <summary>
	/// Private method
	/// Interpolates an influence value utilizing NavMesh_CellGenerator's GrowthValue and the given maxInfluence
	/// </summary>
	/// <param name="maxInfluence">Maximum influence from all neighboring cells</param>
	/// <returns>The newly interpolated influence value for this cell</returns>
	private float CalculateNewInfluence(float maxInfluence)
	{
		return (1.0f - parent.GrowthFactor) * influenceValue + parent.GrowthFactor * maxInfluence;
	}

	/// <summary>
	/// Spawns a new Cell object
	/// A new Dot will be instantiated at the given site with the given rotation
	/// </summary>
	/// <param name="_site">The world position of this cell</param>
	/// <param name="_rotation">The rotation of this cell</param>
	/// <param name="_dotPrefab">The Prefab Dot to be spawned</param>
	/// <param name="_parent">NavMesh_CellGenerator object that this cell belongs to</param>
	public Cell(Vector3 _site, Quaternion _rotation, GameObject _dotPrefab, NavMesh_CellGenerator _parent)
	{
		site = _site;
		influenceValue = 0.0f;
		currentFadeTime = 0.0f;

		parent = _parent;

		dot = UnityEngine.Object.Instantiate(_dotPrefab, site, _rotation);
		dotRenderer = dot.GetComponent<Renderer>();
		dotRenderer.material.SetColor("_Color", Color.white);

		adjacentCells = new List<Pair<Cell, float>>();
	}

	/// <summary>
	/// Given an index value, calls this cell's dot object's PoissonDot.SetParentAndIndex with 
	/// this cell's parent and the given index
	/// </summary>
	/// <param name="index">An index into NavMesh_CellGenerator's cells list</param>
	public void SetIndex(int index)
	{
		dot.GetComponent<PoissonDot>().SetParentAndIndex(parent, index);
	}

	public void MoveSiteVertically(Vector3 verticalOffset)
	{
		site += verticalOffset;
		dot.transform.position = site;
	}

	/// <summary>
	/// Sets this cell's influenceValue to the given value and calls ColorDot
	/// </summary>
	/// <param name="_influenceValue">This cell's new influenceValue</param>
	public void SetInfluenceValue(float _influenceValue)
	{
		influenceValue = _influenceValue;
		ColorDot();
	}

	/// <summary>
	/// Resets currentFadeTime to MaxFadeTimer and sets influenceValue to the given value
	/// </summary>
	/// <param name="_influenceValue">This cell's new influenceValue</param>
	public void SetInfluenceValueAndResetTimer(float _influenceValue)
	{
		influenceValue = _influenceValue;
		currentFadeTime = influenceValue * parent.MaxFadeTimer;
	}

	/// <summary>
	/// For all neighboring cells, if the cell's influenceValue is less than the new
	/// influence value, that cell's SetInfluenceValueAndResetTimer method is called
	/// </summary>
	/// <param name="_influenceValue">influenceValue to pass to neighboring cells</param>
	public void SetNeighborInfluenceValuesAndResetTimer(float _influenceValue)
	{
		foreach (Pair<Cell, float> neighbor in adjacentCells)
		{
			if (neighbor.first.influenceValue < _influenceValue)
				neighbor.first.SetInfluenceValueAndResetTimer(_influenceValue);
		}
	}

	/// <summary>
	/// Sets this cell's influenceValue to this cell's newInfluenceValue and calls ColorDot
	/// </summary>
	public void ApplyNewInfluence()
	{
		influenceValue = newInfluenceValue;
		ColorDot();
	}

	/// <summary>
	/// Given the maximum negative and positive values of all cells, 
	/// this method will divide this cell's newInfluenceValue by the 
	/// respective maximum before calling ApplyNewInfluence
	/// </summary>
	/// <param name="maxNeg">Maximum negative value among all cells (value should be positive)</param>
	/// <param name="maxPos">Maximum positive value among all cells</param>
	public void NormalizeAndApplyNewInfluence(float maxNeg, float maxPos)
	{
		if (newInfluenceValue < 0.0f)
			newInfluenceValue /= maxNeg;
		else
			newInfluenceValue /= maxPos;

		ApplyNewInfluence();
	}

	/// <summary>
	/// Used in propagation, uses a decay + growth formula to calculate this cell's newInfluenceValue
	/// </summary>
	public void CalculateCurrentInfluence()
	{
		float maxInfluence = 0.0f;
		for (int i = 0; i < adjacentCells.Count; ++i)
		{
			float influence = GetInfluenceFromNeighbor(adjacentCells[i]);
			maxInfluence = Math.Abs(maxInfluence) < Math.Abs(influence) ? influence : maxInfluence;
		}

		newInfluenceValue = CalculateNewInfluence(maxInfluence);
	}

	/// <summary>
	/// Returns newInfluenceValue
	/// </summary>
	/// <returns>This cell's newInfluenceValue field, a float value from -1 to 1</returns>
	public float GetNewInfluenceValue()
	{
		return newInfluenceValue;
	}

	/// <summary>
	/// Decrements currentFadeTime and sets influenceValue to (currentFadeTime / MaxFadeTimer) before calling ColorDot
	/// </summary>
	/// <param name="deltaTime">Time since last call</param>
	public void FadeInfluence(float deltaTime)
	{
		currentFadeTime -= deltaTime;
		if (currentFadeTime < 0) currentFadeTime = 0.0f;
		influenceValue = currentFadeTime / parent.MaxFadeTimer;
		ColorDot();
	}

	/// <summary>
	/// Sets this cell's wallInfluence to the given value
	/// </summary>
	/// <param name="_wallInfluence">This cell's new wallInfluence</param>
	public void SetWallInfluence(float _wallInfluence)
	{
		wallInfluence = _wallInfluence;
	}

	/// <summary>
	/// 
	/// </summary>
	public void ApplyWallInfluence()
	{
		influenceValue = wallInfluence;
		ColorDot();
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
		Vector3 verticalOffsetFinal = Vector3.up * verticalOffset;

		line = UnityEngine.Object.Instantiate(linePrefab);
		LineRenderer lineRenderer = line.GetComponent<LineRenderer>();

		int samples = 2 + (int)(Vector3.Distance(start, end) * samplesPerLength);
		lineRenderer.positionCount = samples;
		Vector3[] linePositions = new Vector3[samples];
		for (int i = 0; i < samples; ++i)
		{
			linePositions[i] = Vector3.Lerp(start, end, (float)i / (samples - 1));
			NavMeshHit hit;
			NavMesh.SamplePosition(linePositions[i], out hit, 100, -1);
			linePositions[i] = hit.position + verticalOffsetFinal;
		}
		lineRenderer.SetPositions(linePositions);

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
		for (int i = 0; i < cells.Count; ++i)
		{
			cells[i].SetIndex(i);
		}
	}


	public void RaiseCells(Vector3 verticalOffset)
	{
		foreach (Cell cell in cells)
			cell.MoveSiteVertically(verticalOffset);
	}

	/// <summary>
	/// For every cell contained within this list, calls cell.SetInfluenceValue with 0
	/// </summary>
	public void ResetInfluences()
	{
		foreach (Cell cell in cells)
		{
			cell.SetInfluenceValue(0.0f);
		}
	}

	/// <summary>
	/// Iterates through the cells and calls CalculateCurrentInfluence
	/// Once all cells have finished calculating influence, 
	/// it will iterate through the cells and call ApplyNewInfluence
	/// </summary>
	public void PropagateCells()
	{
		foreach (Cell cell in cells)
		{
			cell.CalculateCurrentInfluence();
		}

		foreach (Cell cell in cells)
		{
			cell.ApplyNewInfluence();
		}
	}

	/// <summary>
	/// Iterates through the cells and calls CalculateCurrentInfluence
	/// Once all cells have finished calculating influence, 
	/// it will iterate through the cells and call ApplyNewInfluence
	/// Then, finds the maximum positve and negative values of all cells and normalizes their influence values
	/// </summary>
	public void NormalizePropogateCells()
	{
		float maxNeg = 0.0f;
		float maxPos = 0.0f;

		foreach (Cell cell in cells)
		{
			cell.CalculateCurrentInfluence();
			float ival = cell.GetNewInfluenceValue();
			if (ival < 0.0f)
				maxNeg = ival < maxNeg ? ival : maxNeg;
			else
				maxPos = ival > maxPos ? ival : maxPos;
		}
		maxNeg *= -1.0f;
		if (maxPos == 0.0f) maxPos = 1.0f;
		if (maxNeg == 0.0f) maxNeg = 1.0f;


		foreach (Cell cell in cells)
		{
			cell.NormalizeAndApplyNewInfluence(maxNeg, maxPos);
		}
	}

	/// <summary>
	/// Iterates through the cells and calls FadeInfluence
	/// </summary>
	/// <param name="deltaTime">Time since last call</param>
	public void FadeInfluenceOnCells(float deltaTime)
	{
		foreach (Cell cell in cells)
		{
			cell.FadeInfluence(deltaTime);
		}
	}

	/// <summary>
	/// Given an array of wall GameObjects, iterates over all cells
	/// and calculates their wallInfluence values
	/// </summary>
	/// <param name="walls">Array of GameObjects with the appropriate wall tag</param>
	public void FindWallInfluence(GameObject[] walls)
	{
		foreach (Cell cell in cells)
		{
			int closestWallIndex = 0;
			float closestWallDist = (cell.site - walls[0].transform.position).sqrMagnitude;
			for(int i = 1; i<walls.Length; ++i)
			{

				float thisWallDist = (cell.site - walls[i].transform.position).sqrMagnitude;
				if (thisWallDist < closestWallDist)
				{
					closestWallDist = thisWallDist;
					closestWallIndex = i;
				}
			}
			closestWallDist = (float)Math.Sqrt((double)closestWallDist);
			float iVal = 1.0f / (closestWallDist * closestWallDist);
			Mathf.Clamp(iVal, 0.0f, 1.0f);
			Debug.Log("Closest Dist: " + closestWallDist + " :: iVal: " + iVal);
			cell.SetWallInfluence(1.0f);
		}
	}

	/// <summary>
	/// 
	/// </summary>
	public void ApplyWallInfluences()
	{
		foreach (Cell cell in cells)
		{
			cell.ApplyWallInfluence();
		}
	}

	/// <summary>
	/// 
	/// </summary>
	public void FindVisibilityInfluence()
	{
		for (int i = 0; i < cells.Count; ++i)
		{
			for (int j = 0; j < cells.Count; ++j)
			{
				if (i == j) continue;

			}
		}
	}
}

public class NavMesh_CellGenerator : MonoBehaviour
{
	public GameObject InfluenceMapModeText;
	public GameObject NavMeshVertPrefab;                        // Prefab for the object used to display NavMesh vertices
    public GameObject PoissonDotPrefab;                         // Prefab for the object used to display Voronoi site locations
    public GameObject LinePrefab;                               // Prefab for the LineRenderer used to display lines and edges
    public GameObject CellConnectionPrefab;                     // Prefab for the LineRenderer used to display edges between CellList
    public List<Triangle> NavMeshTris = new List<Triangle>();   // All Triangles on the NavMesh in the scene
    public float PoissonTolerance = 0.5f;                       // What is the minimum distance apart each poisson point should be
    public float LineWidth = 0.01f;                             // How wide to draw all the lines in the scene
    public SortedCellList CellList;                         // List of all existing CellList within the scene -- will be used for the influence map calculations
	public InfluenceMode Mode = InfluenceMode.Propagation;		// Current mode for the influence map
	public float DecayFactor = 0.2f;							// Decay factor for propagation formula
	public float GrowthFactor = 0.3f;                           // Growth factor for propagation formula
	public float MaxFadeTimer = 4.0f;							// How many seconds it'll take for a given cell to fade it's influence from 1 to 0, value cannot be 0
	public bool ShouldRenderNeighborLines = false;              // Determines whether or not voronoi neighbor lines should be rendered
	public float VerticalOffset = 0.15f;
	Vector3 verticalOffsetFinal;
	public float StepTime = 0.1f;                               // How long between propagation calculations
	float CurrentTime = 0.0f;                                   // Timer for tracking propagation step time

	/// <summary>
	/// Given a Triangle and a tolerance, the triangle's CellList member will be 
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
            //CellList.AddCell(new Cell(points[points.Count - 1], transform.rotation, PoissonDotPrefab, LinePrefab));
            //CellList.GetRange(points[points.Count - 1], PoissonTolerance * 2);
        }

        // Test each point against the selected final points
        // If it meets the tolerance requirement, add it to pointsFinal
        // else, cull this point
        foreach (Vector3 point in points)
        {
            bool tooClose = false;
            Pair<int, int> bounds = CellList.GetRange(point, PoissonTolerance * 2);
            for (int i = bounds.first; i < bounds.second; ++i)
            {
                if (Vector3.Distance(point, CellList.GetCell(i).site) < tolerance)
                {
                    tooClose = true;
                    break;
                }
            }

            if (!tooClose)
            {
                CellList.AddCell(new Cell(point, transform.rotation, PoissonDotPrefab, this));
                pointsFinal.Add(point);
            }
        }
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
        Pair<int, int> bounds = CellList.GetRange(cell.site, PoissonTolerance * 10);

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
            float dist = Vector3.Distance(cell.site, CellList.GetCell(i).site);
            if (dist < PoissonTolerance || Mathf.Abs(cell.site.y - CellList.GetCell(i).site.y) > PoissonTolerance)
            {
                continue;
            }
            if (dist < PoissonTolerance * 3)
            {
                possibleCells.Add(CellList.GetCell(i));
                if (dist < nearestDist)
                {
                    nearest = possibleCells.Count - 1;
                    nearestDist = dist;
                }
            }
            allCells.Add(CellList.GetCell(i));
        }

        if (possibleCells.Count == 0)
        {
            Debug.Log("No nearby points after distance check!!");
            return;
        }

        //Debug line to the nearest cell
        //GameObject l = Instantiate(LinePrefab);
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
        //lr.widthMultiplier = LineWidth;
        //lr.endWidth = LineWidth / 2;

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

			GameObject line = Instantiate(CellConnectionPrefab);
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
				linePositions[i] = hit.position + verticalOffsetFinal;
            }
            lineRenderer.SetPositions(linePositions);

            lineRenderer.enabled = ShouldRenderNeighborLines;
            lineRenderer.widthMultiplier = LineWidth;
            lineRenderer.endWidth = LineWidth / 4;
            cell.adjacentCells.Add(new Pair<Cell, float>(c, Vector3.Distance(c.site, cell.site))); //Add
        }
    }

	/// <summary>
	/// Given an index into CellList, sets that Cell's influence value
	/// </summary>
	/// <param name="index">Index into CellList List</param>
	/// <param name="influence">Influence value to set this cell to</param>
	public void PropagationFromCell(int index, float influence)
	{
		CellList.GetCell(index).SetInfluenceValue(influence);
	}

	/// <summary>
	/// Given an index into CellList, sets the Cell's influence value to maximum
	/// </summary>
	/// <param name="index">Index into CellList List</param>
	public void InFieldOfView(int index)
	{
		if (Mode == InfluenceMode.FieldOfView)
		{
			CellList.GetCell(index).SetInfluenceValueAndResetTimer(1.0f);
			CellList.GetCell(index).SetNeighborInfluenceValuesAndResetTimer(0.5f);
		}
	}

	// Use this for initialization
	void Start()
    {
		verticalOffsetFinal = Vector3.up * VerticalOffset;
		CellList = new SortedCellList();
		/// Cycling through all child objects, finding OffMeshLinks
		/// Hope is to connect our Cells together via the OffMeshLinks
		//foreach (Transform child in transform)
		//{
		//	OffMeshLink offMeshLink = child.GetComponent<OffMeshLink>();
		//	if (offMeshLink && offMeshLink.enabled)
		//	{
		//		NavMesh.areas
		//		offMeshLink.area
		//	}
		//}

        // Get all NavMesh information
        NavMeshTriangulation tris = NavMesh.CalculateTriangulation();

        // Spawn vertex prefabs at the vertices of the NavMesh
        foreach (Vector3 vert in tris.vertices)
        {
            Instantiate(NavMeshVertPrefab, vert, transform.rotation);
        }

        // Using the indicies of the NavMesh, create Triangles and call PoissonDiscDistribution with it to create
        // Voronoi regions
        for (int index = 0; index < tris.indices.Length; index += 3)
        {
            Triangle tri = new Triangle(
                tris.vertices[tris.indices[index + 0]],
                tris.vertices[tris.indices[index + 1]],
                tris.vertices[tris.indices[index + 2]],
                PoissonDotPrefab, LinePrefab, LineWidth
            );
            NavMeshTris.Add(tri);
            PoissonDiscDistribution(tri, PoissonTolerance);
        }

        for (int i = 0; i < CellList.GetCount(); ++i)
        {
            GetNeighbors(CellList.GetCell(i));
        }
		CellList.RaiseCells(verticalOffsetFinal);
		CellList.AssignIndices();
		CellList.FindWallInfluence(GameObject.FindGameObjectsWithTag("Wall"));
		InfluenceMapModeText.GetComponent<ModeUI>().ModeChange(Mode);

		if (Mode == InfluenceMode.OpennessClosestWall)
			CellList.ApplyWallInfluences();
	}

    // Update is called once per frame
    void Update()
    {
		CurrentTime += Time.deltaTime;
		if (CurrentTime >= StepTime)
		{
			switch (Mode)
			{
				case InfluenceMode.Propagation:
				{
					CellList.PropagateCells();
					break;
				}
				case InfluenceMode.NormalizedPropagation:
				{
					CellList.PropagateCells();
					CellList.NormalizePropogateCells();
					break;
				}
				case InfluenceMode.FieldOfView:
				{
					CellList.FadeInfluenceOnCells(StepTime);
					break;
				}
			}

			CurrentTime = 0.0f;
		}


		if (Input.GetKeyUp(KeyCode.LeftBracket))
		{
			CellList.ResetInfluences();
			Mode = --Mode >= 0 ? Mode : InfluenceMode.NumValues - 1;
			InfluenceMapModeText.GetComponent<ModeUI>().ModeChange(Mode);

			if (Mode == InfluenceMode.OpennessClosestWall)
				CellList.ApplyWallInfluences();
		}
		if (Input.GetKeyUp(KeyCode.RightBracket))
		{
			CellList.ResetInfluences();
			Mode = ++Mode < InfluenceMode.NumValues ? Mode : 0;
			InfluenceMapModeText.GetComponent<ModeUI>().ModeChange(Mode);

			if (Mode == InfluenceMode.OpennessClosestWall)
				CellList.ApplyWallInfluences();
		}
	}
}
