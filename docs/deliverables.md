EC API

STAGE 0: PURPOSE (WHAT DOES IT DO)
* RESTful API
* 

STAGE 1. DIVISION OF RESPONSIBILITY & PLANNING (WHAT DOES IT LOOK LIKE)
Separate projects:
1. FRONT END (API calls, results display, etc.) (CHRIS)
	* What should be presented to user?
	* How should it be presented to user?
	* What set of parameters should be able to be manipulated?
2. BACK END API (Methods, functions, routes, URLS, etc.) (JACKSON)
	* How can raw data be made accessible to end users?
	* How can raw data be manipulated to present new insights to users?
	* How can user inputs be translated into DB queries?
	* How can data be best accessed by third parties (non GEMSEC users)
3. BACK END DB (Schema, storage, ORM)
	* What, exactly, are we using? And the drawbacks/benefits therein
	* What is the exact schema of the database for different sets of data?

STAGE 2. (HOW IS IT CREATED)
1. DOCUMENTATION, DOCUMENTATION, DOCUMENTATION
	* Every specified route should have proper documentation -- largely handled by FastAPI, but make sure to
	  create docstrings in each method adding more detail to usage and limitations

DELIVERABLES
1. Design document (URP)
2. Prototype (simple page to query, fetch, and display data contained in DB) (URP)
3. Schema for database (URP)
4. Service for uploading, storing, viewing data (URP)
5. UW NET ID Authentication (URP), Additional methods of authentication (Spring Quarter).
6. Front End (have a substantially completed version by URP)

7. iPad OS 

friday afternoons are when metadata meetings
tuesday and thursday
Jackson and Chris meeting on Mondays


