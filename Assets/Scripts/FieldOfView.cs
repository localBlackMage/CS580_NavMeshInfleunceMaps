using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FieldOfView : MonoBehaviour {

	// Use this for initialization
	void Start () {
		
	}
	
	// Update is called once per frame
	void Update () {
		
	}

	void OnTriggerStay(Collider other)
	{
		PoissonDot pDot = other.GetComponent<PoissonDot>();
		if (pDot) {
			//Debug.DrawLine(transform.parent.position, other.transform.position, Color.green);
			pDot.InFieldOfView();
		}
	}
}
