using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PoissonDot : MonoBehaviour {

	NavMesh_CellGenerator parent;
    MeshRenderer meshRenderer;
    int index;
    bool rendered = true;

	public void SetParentAndIndex(NavMesh_CellGenerator _parent, int _index)
	{
		parent = _parent;
		index = _index;
	}

	public void InFieldOfView()
	{
		parent.InFieldOfView(index);
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

    private void Start()
    {
        meshRenderer = GetComponent<MeshRenderer>();
    }

    private void Update()
    {
        if (Input.GetKeyUp(KeyCode.M))
        {
            rendered = !rendered;
            meshRenderer.enabled = rendered;
        }
    }
}
