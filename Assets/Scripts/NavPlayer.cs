using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NavPlayer : MonoBehaviour {

	public float cameraSensitivity = 90;
	public float normalMoveSpeed = 10;

	private float rotationX = 0.0f;
	private float rotationZ = 0.0f;

	private bool locked = true;

	private Transform head;

	// Use this for initialization
	void Start ()
	{
		head = transform.GetChild(0);
	}
	
	// Update is called once per frame
	void Update ()
	{
		if (Input.GetKeyDown(KeyCode.Tab))
		{
			locked = !locked;
		}

		if (!locked)
		{
			rotationX += Input.GetAxis("Mouse X") * cameraSensitivity * Time.deltaTime;
			rotationZ += Input.GetAxis("Mouse Y") * cameraSensitivity * Time.deltaTime;
			rotationZ = Mathf.Clamp(rotationZ, -90, 90);

			head.localRotation = Quaternion.AngleAxis(rotationX, Vector3.up);
			head.localRotation *= Quaternion.AngleAxis(0, Vector3.left);
			head.localRotation *= Quaternion.AngleAxis(rotationZ, Vector3.forward);

			transform.position += transform.forward * normalMoveSpeed * Input.GetAxis("Vertical") * Time.deltaTime;
			transform.position += transform.right * normalMoveSpeed * Input.GetAxis("Horizontal") * Time.deltaTime;
		}
	}
}
