/**
 * Copyright (C) 2011,2012 Gordon Fraser, Andrea Arcuri and EvoSuite
 * contributors
 * 
 * This file is part of EvoSuite.
 * 
 * EvoSuite is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 * 
 * EvoSuite is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Public License for more details.
 * 
 * You should have received a copy of the GNU Public License along with
 * EvoSuite. If not, see <http://www.gnu.org/licenses/>.
 */
package org.evosuite.testcase;

import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.ClassUtils;
import org.evosuite.Properties;
import org.evosuite.utils.GenericClass;
import org.evosuite.utils.GenericConstructor;
import org.evosuite.utils.Randomness;
import org.objectweb.asm.Label;
import org.objectweb.asm.Opcodes;
import org.objectweb.asm.Type;
import org.objectweb.asm.commons.GeneratorAdapter;
import org.objectweb.asm.commons.Method;

import edu.uta.cse.dsc.VM;

/**
 * This statement represents a constructor call
 * 
 * @author Gordon Fraser
 */
public class ConstructorStatement extends AbstractStatement {

	private static final long serialVersionUID = -3035570485633271957L;

	private GenericConstructor constructor;

	public List<VariableReference> parameters;

	private static List<String> primitiveClasses = Arrays.asList("char", "int", "short",
	                                                             "long", "boolean",
	                                                             "float", "double",
	                                                             "byte");

	/**
	 * <p>
	 * Constructor for ConstructorStatement.
	 * </p>
	 * 
	 * @param tc
	 *            a {@link org.evosuite.testcase.TestCase} object.
	 * @param constructor
	 *            a {@link java.lang.reflect.Constructor} object.
	 * @param type
	 *            a {@link java.lang.reflect.Type} object.
	 * @param parameters
	 *            a {@link java.util.List} object.
	 */
	public ConstructorStatement(TestCase tc, GenericConstructor constructor,
	        List<VariableReference> parameters) {
		super(tc, new VariableReferenceImpl(tc, constructor.getOwnerClass()));
		this.constructor = constructor;
		// this.return_type = constructor.getDeclaringClass();
		this.parameters = parameters;
	}

	/**
	 * This constructor allows you to use an already existing VariableReference
	 * as retvar. This should only be done, iff an old statement is replaced
	 * with this statement. And already existing objects should in the future
	 * reference this object.
	 * 
	 * @param tc
	 *            a {@link org.evosuite.testcase.TestCase} object.
	 * @param constructor
	 *            a {@link java.lang.reflect.Constructor} object.
	 * @param retvar
	 *            a {@link org.evosuite.testcase.VariableReference} object.
	 * @param parameters
	 *            a {@link java.util.List} object.
	 */
	public ConstructorStatement(TestCase tc, GenericConstructor constructor,
	        VariableReference retvar, List<VariableReference> parameters) {
		super(tc, retvar);
		assert (tc.size() > retvar.getStPosition()); //as an old statement should be replaced by this statement
		this.constructor = constructor;
		// this.return_type = constructor.getDeclaringClass();
		this.parameters = parameters;
	}

	/**
	 * <p>
	 * Constructor for ConstructorStatement.
	 * </p>
	 * 
	 * @param tc
	 *            a {@link org.evosuite.testcase.TestCase} object.
	 * @param constructor
	 *            a {@link java.lang.reflect.Constructor} object.
	 * @param retvar
	 *            a {@link org.evosuite.testcase.VariableReference} object.
	 * @param parameters
	 *            a {@link java.util.List} object.
	 * @param check
	 *            a boolean.
	 */
	protected ConstructorStatement(TestCase tc, GenericConstructor constructor,
	        VariableReference retvar, List<VariableReference> parameters, boolean check) {
		super(tc, retvar);
		assert check == false;
		this.constructor = constructor;
		this.parameters = parameters;
	}

	/**
	 * <p>
	 * Getter for the field <code>constructor</code>.
	 * </p>
	 * 
	 * @return a {@link java.lang.reflect.Constructor} object.
	 */
	public GenericConstructor getConstructor() {
		return constructor;
	}

	/**
	 * <p>
	 * Setter for the field <code>constructor</code>.
	 * </p>
	 * 
	 * @param constructor
	 *            a {@link java.lang.reflect.Constructor} object.
	 */
	public void setConstructor(GenericConstructor constructor) {
		this.constructor = constructor;
		retval.setType(constructor.getReturnType());
	}

	/**
	 * <p>
	 * getReturnType
	 * </p>
	 * 
	 * @param clazz
	 *            a {@link java.lang.Class} object.
	 * @return a {@link java.lang.String} object.
	 */
	public static String getReturnType(Class<?> clazz) {
		String retVal = ClassUtils.getShortClassName(clazz);
		if (primitiveClasses.contains(retVal))
			return clazz.getSimpleName();

		return retVal;
	}

	// TODO: Handle inner classes (need instance parameter for newInstance)
	/** {@inheritDoc} */
	@Override
	public Throwable execute(final Scope scope, PrintStream out)
	        throws InvocationTargetException, IllegalArgumentException,
	        InstantiationException, IllegalAccessException {
		//PrintStream old_out = System.out;
		//PrintStream old_err = System.err;
		//System.setOut(out);
		//System.setErr(out);

		logger.trace("Executing constructor " + constructor.toString());
		final Object[] inputs = new Object[parameters.size()];

		try {
			return super.exceptionHandler(new Executer() {

				@Override
				public void execute() throws InvocationTargetException,
				        IllegalArgumentException, IllegalAccessException,
				        InstantiationException, CodeUnderTestException {

					java.lang.reflect.Type[] parameterTypes = constructor.getParameterTypes();
					for (int i = 0; i < parameters.size(); i++) {
						VariableReference parameterVar = parameters.get(i);
						if (!parameterVar.isAssignableTo(parameterTypes[i])) {
							throw new CodeUnderTestException(
							        new UncompilableCodeException());
						}
						try {
							inputs[i] = parameterVar.getObject(scope);
						} catch (CodeUnderTestException e) {
							throw e;
							//throw new CodeUnderTestException(e.getCause());
							// throw CodeUnderTestException.throwException(e.getCause());
						} catch (Throwable e) {
							//FIXME: this does not seem to propagate to client root. Is this normal behavior?
							logger.error("Class " + Properties.TARGET_CLASS
							        + ". Error encountered: " + e);
							assert (false);
							throw new EvosuiteError(e);
						}
					}

					// If this is a non-static member class, the first parameter must not be null
					if (constructor.getConstructor().getDeclaringClass().isMemberClass()
					        && !Modifier.isStatic(constructor.getConstructor().getDeclaringClass().getModifiers())) {
						if (inputs[0] == null) {
							// throw new NullPointerException();
							throw new CodeUnderTestException(new NullPointerException());
						}
					}

					Object ret = constructor.getConstructor().newInstance(inputs);

					try {
						// assert(retval.getVariableClass().isAssignableFrom(ret.getClass())) :"we want an " + retval.getVariableClass() + " but got an " + ret.getClass();
						retval.setObject(scope, ret);
					} catch (CodeUnderTestException e) {
						throw e;
						// throw CodeUnderTestException.throwException(e);
					} catch (Throwable e) {
						throw new EvosuiteError(e);
					}
				}

				@Override
				public Set<Class<? extends Throwable>> throwableExceptions() {
					Set<Class<? extends Throwable>> t = new HashSet<Class<? extends Throwable>>();
					t.add(InvocationTargetException.class);
					return t;
				}
			});

		} catch (InvocationTargetException e) {
			VM.setIgnoreCallBack(true);
			exceptionThrown = e.getCause();
			if (logger.isDebugEnabled()) {
				try {
					logger.debug("Exception thrown in constructor: " + e.getCause());
				}
				//this can happen if SUT throws exception on toString
				catch (Exception ex) {
					logger.debug("Exception thrown in constructor and SUT gives issue when calling e.getCause()",
					             ex);
				}
			}
		}
		return exceptionThrown;
	}

	/** {@inheritDoc} */
	@Override
	public StatementInterface copy(TestCase newTestCase, int offset) {
		ArrayList<VariableReference> new_params = new ArrayList<VariableReference>();
		for (VariableReference r : parameters) {
			new_params.add(r.copy(newTestCase, offset));
		}

		AbstractStatement copy = new ConstructorStatement(newTestCase, constructor,
		        new_params);
		// copy.assertions = copyAssertions(newTestCase, offset);

		return copy;
	}

	/** {@inheritDoc} */
	@Override
	public Set<VariableReference> getVariableReferences() {
		Set<VariableReference> references = new LinkedHashSet<VariableReference>();
		references.add(retval);
		references.addAll(parameters);
		for (VariableReference param : parameters) {
			if (param.getAdditionalVariableReference() != null)
				references.add(param.getAdditionalVariableReference());
		}
		references.addAll(getAssertionReferences());

		return references;
	}

	/* (non-Javadoc)
	 * @see org.evosuite.testcase.StatementInterface#replace(org.evosuite.testcase.VariableReference, org.evosuite.testcase.VariableReference)
	 */
	/** {@inheritDoc} */
	@Override
	public void replace(VariableReference var1, VariableReference var2) {

		if (retval.equals(var1))
			retval = var2;

		for (int i = 0; i < parameters.size(); i++) {

			if (parameters.get(i).equals(var1))
				parameters.set(i, var2);
			else
				parameters.get(i).replaceAdditionalVariableReference(var1, var2);
		}
	}

	/**
	 * <p>
	 * getParameterReferences
	 * </p>
	 * 
	 * @return a {@link java.util.List} object.
	 */
	public List<VariableReference> getParameterReferences() {
		return parameters;
	}

	/* (non-Javadoc)
	 * @see org.evosuite.testcase.StatementInterface#getNumParameters()
	 */
	@Override
	public int getNumParameters() {
		return parameters.size();
	}

	/**
	 * <p>
	 * replaceParameterReference
	 * </p>
	 * 
	 * @param var
	 *            a {@link org.evosuite.testcase.VariableReference} object.
	 * @param numParameter
	 *            a int.
	 */
	public void replaceParameterReference(VariableReference var, int numParameter) {
		assert (numParameter >= 0);
		assert (numParameter < parameters.size());
		parameters.set(numParameter, var);
	}

	/** {@inheritDoc} */
	@Override
	public boolean equals(Object s) {
		if (this == s)
			return true;
		if (s == null)
			return false;
		if (getClass() != s.getClass())
			return false;

		ConstructorStatement ms = (ConstructorStatement) s;
		if (ms.parameters.size() != parameters.size())
			return false;

		if (!this.constructor.equals(ms.constructor))
			return false;

		for (int i = 0; i < parameters.size(); i++) {
			if (!parameters.get(i).equals(ms.parameters.get(i)))
				return false;
		}

		return retval.equals(ms.retval);
	}

	/** {@inheritDoc} */
	@Override
	public int hashCode() {
		final int prime = 41;
		int result = 1;
		result = prime * result + ((constructor == null) ? 0 : constructor.hashCode());
		result = prime * result + ((parameters == null) ? 0 : parameters.hashCode());
		return result;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.evosuite.testcase.Statement#getBytecode(org.objectweb.
	 * asm.commons.GeneratorAdapter)
	 */
	/** {@inheritDoc} */
	@Override
	public void getBytecode(GeneratorAdapter mg, Map<Integer, Integer> locals,
	        Throwable exception) {
		logger.debug("Invoking constructor");
		Label start = mg.newLabel();
		Label end = mg.newLabel();

		// if(exception != null)
		mg.mark(start);

		mg.newInstance(Type.getType(retval.getVariableClass()));
		mg.dup();
		int num = 0;
		for (VariableReference parameter : parameters) {
			parameter.loadBytecode(mg, locals);
			if (constructor.getConstructor().getParameterTypes()[num].isPrimitive()) {
				if (parameter.getGenericClass().isWrapperType()) {
					mg.unbox(Type.getType(parameter.getGenericClass().getUnboxedType()));
				} else if (!parameter.getGenericClass().isPrimitive()) {
					Class<?> parameterClass = new GenericClass(
					        constructor.getParameterTypes()[num]).getBoxedType();
					Type parameterType = Type.getType(parameterClass);
					mg.checkCast(parameterType);
					mg.unbox(Type.getType(constructor.getConstructor().getParameterTypes()[num]));
				}

				if (!constructor.getParameterTypes()[num].equals(parameter.getVariableClass())) {
					logger.debug("Types don't match - casting "
					        + parameter.getVariableClass().getName()
					        + " to "
					        + constructor.getConstructor().getParameterTypes()[num].getName());
					mg.cast(Type.getType(parameter.getVariableClass()),
					        Type.getType(constructor.getConstructor().getParameterTypes()[num]));
				}
			} else if (parameter.getVariableClass().isPrimitive()) {
				mg.box(Type.getType(parameter.getVariableClass()));
			}
			num++;
		}
		mg.invokeConstructor(Type.getType(retval.getVariableClass()),
		                     Method.getMethod(constructor.getConstructor()));
		logger.debug("Storing result");
		retval.storeBytecode(mg, locals);

		// if(exception != null) {
		mg.mark(end);
		Label l = mg.newLabel();
		mg.goTo(l);
		// mg.catchException(start, end,
		// Type.getType(getExceptionClass(exception)));
		mg.catchException(start, end, Type.getType(Throwable.class));
		mg.pop(); // Pop exception from stack
		if (!retval.isVoid()) {
			Class<?> clazz = retval.getVariableClass();
			if (clazz.equals(boolean.class))
				mg.push(false);
			else if (clazz.equals(char.class))
				mg.push(0);
			else if (clazz.equals(int.class))
				mg.push(0);
			else if (clazz.equals(short.class))
				mg.push(0);
			else if (clazz.equals(long.class))
				mg.push(0L);
			else if (clazz.equals(float.class))
				mg.push(0.0F);
			else if (clazz.equals(double.class))
				mg.push(0.0);
			else if (clazz.equals(byte.class))
				mg.push(0);
			else if (clazz.equals(String.class))
				mg.push("");
			else
				mg.visitInsn(Opcodes.ACONST_NULL);

			retval.storeBytecode(mg, locals);
		}
		mg.mark(l);
		// }

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.evosuite.testcase.Statement#getDeclaredExceptions()
	 */
	/** {@inheritDoc} */
	@Override
	public Set<Class<?>> getDeclaredExceptions() {
		Set<Class<?>> ex = super.getDeclaredExceptions();
		ex.addAll(Arrays.asList(constructor.getConstructor().getExceptionTypes()));
		return ex;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.evosuite.testcase.Statement#getUniqueVariableReferences()
	 */
	/** {@inheritDoc} */
	@Override
	public List<VariableReference> getUniqueVariableReferences() {
		List<VariableReference> references = new ArrayList<VariableReference>();
		references.add(retval);
		references.addAll(parameters);
		for (VariableReference param : parameters) {
			if (param instanceof ArrayIndex)
				references.add(((ArrayIndex) param).getArray());
		}
		return references;

	}

	/**
	 * Go through parameters of constructor call and apply local search
	 * 
	 * @param test
	 * @param statement
	 * @param objective
	 */
	/* (non-Javadoc)
	 * @see org.evosuite.testcase.AbstractStatement#mutate(org.evosuite.testcase.TestCase, org.evosuite.testcase.TestFactory)
	 */
	@Override
	public boolean mutate(TestCase test, TestFactory factory) {

		if (Randomness.nextDouble() >= Properties.P_CHANGE_PARAMETER)
			return false;

		List<VariableReference> parameters = getParameterReferences();
		if (parameters.isEmpty())
			return false;

		int numParameter = Randomness.nextInt(parameters.size());
		VariableReference parameter = parameters.get(numParameter);

		List<VariableReference> objects = test.getObjects(parameter.getType(),
		                                                  getPosition());
		objects.remove(parameter);
		objects.remove(getReturnValue());

		NullStatement nullStatement = new NullStatement(test, parameter.getType());
		StatementInterface primitiveCopy = null;

		if (!parameter.isPrimitive())
			objects.add(nullStatement.getReturnValue());
		else {
			StatementInterface originalStatement = test.getStatement(parameter.getStPosition());
			if (originalStatement instanceof PrimitiveStatement<?>) {
				PrimitiveStatement<?> copy = (PrimitiveStatement<?>) originalStatement.clone(test);
				copy.delta();
				objects.add(copy.getReturnValue());
				primitiveCopy = copy;
			}
		}

		if (objects.isEmpty())
			return false;

		VariableReference replacement = Randomness.choice(objects);
		if (replacement == nullStatement.getReturnValue()) {
			test.addStatement(nullStatement, getPosition());
		} else if (primitiveCopy != null && replacement == primitiveCopy.getReturnValue()) {
			test.addStatement(primitiveCopy, getPosition());
		}

		replaceParameterReference(replacement, numParameter);

		return true;
	}

	/* (non-Javadoc)
	 * @see org.evosuite.testcase.StatementInterface#isValid()
	 */
	/** {@inheritDoc} */
	@Override
	public boolean isValid() {
		assert (super.isValid());
		for (VariableReference v : parameters) {
			v.getStPosition();
		}
		return true;
	}

	/** {@inheritDoc} */
	@Override
	public boolean same(StatementInterface s) {
		if (this == s)
			return true;
		if (s == null)
			return false;
		if (getClass() != s.getClass())
			return false;

		ConstructorStatement ms = (ConstructorStatement) s;
		if (ms.parameters.size() != parameters.size())
			return false;

		if (!this.constructor.equals(ms.constructor))
			return false;

		for (int i = 0; i < parameters.size(); i++) {
			if (!parameters.get(i).same(ms.parameters.get(i)))
				return false;
		}

		return retval.same(ms.retval);
	}

	/** {@inheritDoc} */
	@Override
	public GenericConstructor getAccessibleObject() {
		return constructor;
	}

	/** {@inheritDoc} */
	@Override
	public boolean isAssignmentStatement() {
		return false;
	}

	/* (non-Javadoc)
	 * @see org.evosuite.testcase.StatementInterface#changeClassLoader(java.lang.ClassLoader)
	 */
	/** {@inheritDoc} */
	@Override
	public void changeClassLoader(ClassLoader loader) {
		constructor.changeClassLoader(loader);
		super.changeClassLoader(loader);
	}
}
