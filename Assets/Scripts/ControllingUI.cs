using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public enum ControlMode : int
{
	Camera = 0,
	Player
}

public class ControllingUI : MonoBehaviour {

	Text textRenderer;

	public string[] ControlStrings = {
		"Camera",
		"Player"
	};
	public string BaseText = "Controlling: ";	// Base string for this text component

	public void ControlChange(ControlMode mode)
	{
		if (!textRenderer)
			textRenderer = GetComponent<Text>();
		textRenderer.text = BaseText + ControlStrings[(int)mode];
	}

	// Use this for initialization
	void Start () {
		textRenderer = GetComponent<Text>();
	}
	
	// Update is called once per frame
	void Update () {
		
	}
}
