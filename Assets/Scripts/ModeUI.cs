using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ModeUI : MonoBehaviour {
	Text textRenderer;

	public string[] InfluenceModeStrings = {
		"Propagation",
		"Normalized Propagation",
		"Field of View",
		"Visible to Cell",
		"Openness/Closest Wall"
	};
	public GameObject TerrainParent;					// Nav Mesh terrain in world
	public string BaseText = "Influence Map Mode: ";	// Base string for this text component

	public void ModeChange(InfluenceMode mode)
	{
		if (!textRenderer)
			textRenderer = GetComponent<Text>();
		textRenderer.text = BaseText + InfluenceModeStrings[(int)mode];
	}

	// Use this for initialization
	void Start () {
		textRenderer = GetComponent<Text>();
		textRenderer.text = BaseText + InfluenceModeStrings[
			(int)TerrainParent.GetComponent<NavMesh_CellGenerator>().Mode
		];
	}
	
	// Update is called once per frame
	void Update () {
		
	}
}
