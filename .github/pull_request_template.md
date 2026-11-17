*Replace with one-line Summary*

Sci/Tech Reviewer: <!-- SR id, filled by author when ready for review -->  
Code Reviewer: <!-- CR id, filled by system team -->
<!-- IAO/IAM: include-if-required -->

---

## Description <!-- To be completed by the developer -->
<!-- Provide a brief description of the changes in this PR, including any notes
     useful for reviewers -->

### :link: Related Issues/PRs <!-- delete this section for standalone PRs -->
<!-- List any related PRs here
- issue/prs from the same repository can be added by #Number
- issues/prs from other repositories can be included as Org/Repo#Number
-->

- closes #issue-number (auto-closes the issue)
- fixes #issue-number (auto-closes the issue)
- is related to #issue-number, #pr-number
- linked #pr-number
- blocked-by #pr-number
- blocks #pr-number

### :memo: Contributor Checklist

- [ ] I have performed a self-review of my own code.
- [ ] All commits in this PR are signed and show as verified.
- [ ] My code follows the project's [style guidelines](https://metoffice.github.io/lfric_core/how_to_contribute/index.html#how-to-contribute-index).
- [ ] Comments have been included that aid understanding and enhance the readability of the code.
- [ ] Where appropriate I have updated documentation related to this change and confirmed that it builds correctly.
- [ ] My changes generate no new warnings.

### :test_tube: Testing

- [ ] I have tested this change locally, using the LFRic Apps rose-stem suite.
- [ ] If any tests fail (rose-stem or CI) the reason is understood and acceptable (e.g., KGO changes).
- [ ] I have added tests to cover new functionality as appropriate (e.g. system tests, unit tests, etc.).
- [ ] Any new tests have been assigned an appropriate amount of compute resource and have been allocated to an appropriate testing group (i.e., the developer tests are for jobs which use a small amount of compute resource and complete in a matter of minutes).

<!-- Describe other testing performed (if applicable) -->

#### trac.log

<!-- Paste your trac.log from testing output here -->

### :shield: Security Considerations

- [ ] I have reviewed my changes for potential security issues.
- [ ] Sensitive data is properly handled (if applicable).
- [ ] Authentication and authorisation are properly implemented (if applicable).

### :zap: Performance Impact

- [ ] Performance of the code has been considered and, if applicable, suitable performance measurements have been conducted.  

### :robot: Generative AI Declaration

**Any Generative AI tools (e.g., GitHub Copilot, Calude Code, ChatGPT) to assist in writing, refactoring, or optimising code for this PR?**

- [ ] **No.** This code was written entirely manually.
- [ ] **Yes.** (If yes, please complete the verification items below).

<details>
<summary><b>Click here to expand the AI Usage Verification Checklist</b></summary>

#### AI Usage Verification

*Only check the boxes below if they are true. If any box cannot be checked, your PR will be automatically rejected.*

- [ ] **Approved Tool:** The AI tool used is an approved enterprise/institutional tier model that provides IPR indemnity. An exception may be made with IAO/IAM approval, subject to appropriate scrutiny of the generated code.  
`Assisted-by: [Institution Name] [Tool / Agent Context] ([Model Version]).`  
See [AI attribution requirements](https://metoffice.github.io/simulation-systems/FurtherDetails/ai.html#attribution-requirements).
- [ ] **No Free/Personal Tiers:** I confirm that **no** free or personal tier AI tools (which lack legal protection and IPR indemnity) were used.
- [ ] **In-File Attribution:** I have added the required `! Some content in this file was generated...` comment block following the [Simulation Systems AI policy](https://metoffice.github.io/simulation-systems/FurtherDetails/ai.html).
Commit-message attribution is encouraged for traceability, but is optional. It does not replace the required attribution in the PR description and source file headers.
- [ ] **Code Responsibility:** I have personally reviewed, tested, and understood all AI-assisted sections, and I take full responsibility for its logic, security, and IPR compliance.

</details>

### :white_check_mark: Specialist Approval

- [ ] If this PR includes optimisation work or is expected to affect performance, I have obtained approval from `@MetOffice/HPC-opt` team.
- [ ] If this PR includes PSyclone-related changes (e.g. PSyKAl-lite, Kernel interface, optimisation scripts, LFRic data structure code), I have obtained approval from `@christophermaynard`.
- [ ] If this PR includes XIOS-related changes, I have obtained approval from `@MetOffice/ToolsCollabDev` team.

## :microscope: Sci/Tech Review <!-- To be completed by the Sci/Tech Reviewer -->
<!-- May be skipped for trivial tickets -->

- [ ] I have reviewed the scientific approach and its implementation described in the PR and consider them appropriate.
- [ ] I have identified whether scientific documentation or papers need updating.

(_Please alert the code reviewer via a tag when you have approved the SR_)

## :computer: Code Review <!-- To be completed by the Code Reviewer -->

- [ ] All dependencies have been resolved.
- [ ] All relevant CODEOWNERS have approved the PR, where applicable.
- [ ] Related Issues have been properly linked and addressed.
- [ ] Code quality standards have been met.
- [ ] Tests are adequate and have passed.
- [ ] Documentation is complete and accurate.
- [ ] Security considerations have been addressed.
- [ ] Performance impact is acceptable.
- [ ] PR complies with the [Simulation Systems AI policy](https://metoffice.github.io/simulation-systems/FurtherDetails/ai.html).
