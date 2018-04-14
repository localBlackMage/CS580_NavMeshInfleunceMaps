using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PoissonDot : MonoBehaviour {

	NavMesh_CellGenerator parent;
	int index;

	public void SetParentAndIndex(NavMesh_CellGenerator _parent, int _index)
	{
		parent = _parent;
		index = _index;
	}

	void OnMouseOver()
	{
		if (Input.GetMouseButtonUp(0))
		{
			parent.PropagationFromCell(index, 1.0f);
		}
		else if (Input.GetMouseButtonUp(1))
		{
			parent.PropagationFromCell(index, -1.0f);
		}
	}
}
